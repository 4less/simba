import os
from typing import List, Tuple, Any

import numpy as np
import re
from src.distributions import *
from src.genome_resource import Rank, Genome
from loguru import logger

import random

from src.simulator_wrappers import Simulator
from src.utils import make_executable, longest_common_prefix, keep_only_folder, mkdir_if_not_exists
from src.genome_resource import Genome

def replace_extension(input_string, replacement_extension=''):
    # Define a regular expression pattern for ".fna.gz" or ".fna" at the end of a string
    pattern = re.compile(r'\.fna(\.gz)?$')

    # Use re.sub() to replace the matched pattern with your desired replacement
    result = re.sub(pattern, replacement_extension, input_string)

    return result


def scale_counts(counts, min_vcov, max_vcov):
    min_c = min(counts)
    counts = [c - min_c for c in counts]
    max_c = max(counts)
    
    factor: float = (max_vcov - min_vcov) / max_c
    
    return [(factor * count) + min_vcov for count in counts]


def batch_list(input_list, batch_size: int):
    return [list(range(i, min(i + batch_size, len(input_list)))) for i in range(0, len(input_list), batch_size)]

class Generator:
    @staticmethod
    def generate_vertical_coverage(genome_to_group, min_vcov, max_vcov, method="pareto"):
        """Generates vertical coverages for genomes

        Parameters
        ----------
        genome_to_group : list[]
            groups of genomes to simulate coverages together for
        min_vcov : float
            A flag used to print the columns to the console (default is
            False)
        max_vcov : float
            A flag used to print the columns to the console (default is
            False)
        method : str
            Different methods for generating vertical coverage

        Returns
        -------
        dict
            a list of strings used that are the header columns
        """

        group_list = list(set(genome_to_group.values()))
        random.shuffle(group_list)
        m = len(group_list)
        
        if method == "pareto":
            initial_vcov = pareto_distr(m)
        elif method == "linear":
            initial_vcov = linear_distr(m)
        else:
            initial_vcov = linear_distr(m)
        
        scaled_vcov = scale_counts(initial_vcov, min_vcov, max_vcov)
        


        assert max_vcov == max(scaled_vcov)
        assert min_vcov == min(scaled_vcov)

        group_to_cov = dict(zip(group_list, scaled_vcov))

        assert min_vcov == min(vcov for group, vcov in group_to_cov.items())
        assert max_vcov == max(vcov for group, vcov in group_to_cov.items())

        return group_to_cov
        
    @staticmethod
    def write_simulate_reads_script(output, joint_read_output_folder, simulation_output_folder, genomes, coverages, sample_id, gzip=False):
        output.write('#!/bin/bash\n\n')
        
        # output.write('SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )\n')
        output.write("cd {}\n\n".format(os.getcwd()))

        genomes = list(genomes)
        
        paths = [g.path for g in genomes]
        lcp = longest_common_prefix(paths)
        lcpath = keep_only_folder(lcp)
        lcpath = lcpath + '/' if len(lcpath) > 0 else ""
        
        VARIABLE_STRING="PATH_PREFIX"
        output.write("{}={}\n".format(VARIABLE_STRING, lcpath))
        
        SIMULATION_OUTPUT="READ_OUTPUT_FOLDER"
        output.write("{}={}\n".format(SIMULATION_OUTPUT, simulation_output_folder))

        SIMULATION_OUTPUT_JOINT="JOINT_READ_OUTPUT_FOLDER"
        output.write("{}={}\n".format(SIMULATION_OUTPUT_JOINT, joint_read_output_folder))
        
        
        output.write("mkdir -p {}\n".format(simulation_output_folder))
        output.write("mkdir -p {}\n".format(joint_read_output_folder))

        read_list_fwd = []
        read_list_rev = []

        for genome, coverage in zip(genomes, coverages):
            shortened_genome_path = genome.unzipped_path().replace(lcpath, "") if len(lcpath) > 0 else genome.path
            genome_str = "${{{}}}/{}".format(VARIABLE_STRING, shortened_genome_path).replace('.gz', '')

            shell_command, read_fwd, read_rev = Simulator.get_art_illumina(genome_str, "${{{}}}/{}".format(SIMULATION_OUTPUT, genome.id), True, 150, coverage)

            read_list_fwd.append(read_fwd)
            read_list_rev.append(read_rev)

            output.write(' '.join(map(str, shell_command)))
            output.write('\n')

        read_fwd_out = "${{{}}}/{}_1.fq".format(SIMULATION_OUTPUT_JOINT, sample_id)
        read_rev_out = "${{{}}}/{}_2.fq".format(SIMULATION_OUTPUT_JOINT, sample_id)

        output.write("cat {} > {}\n".format(' '.join(read_list_fwd), read_fwd_out))
        output.write("cat {} > {}\n".format(' '.join(read_list_rev), read_rev_out))

        if gzip:
            output.write("gzip {}\n".format(' '.join(read_list_fwd)))
            output.write("gzip {}\n".format(' '.join(read_list_rev)))
            output.write("gzip {}\n".format(read_fwd_out))
            output.write("gzip {}\n".format(read_rev_out))

        output.write('\n')


    @staticmethod
    def write_conspecific_profile(output, genomes: list[Genome], coverages, delimiter='\t'):
        species = set([g.r_species for g in genomes])
        genome_to_coverage = { g.id: vcov for g, vcov in zip(genomes, coverages) }

        total_coverage = sum(coverages)
        for sp in species:
            sp_genomes: list[tuple[Genome, float]] = [(g, genome_to_coverage[g.id]) for g in genomes if g.r_species == sp]
            sp_total_coverage = sum(vcov for g, vcov in sp_genomes)
            sp_abundance = sp_total_coverage / total_coverage

            sorted(sp_genomes, key=lambda pair: pair[1], reverse=True)

            for genome, coverage in sp_genomes:
                output.write('\t'.join(map(str, [
                    genome.id,
                    genome.r_species,
                    sp_total_coverage,
                    sp_abundance,
                    coverage,
                    coverage / sp_total_coverage,
                    genome.path
                ])))
                output.write('\n')

    @staticmethod
    def write_profile(output, genomes: set[Genome], coverages, delimiter='\t'):
        total_coverage = sum(coverages)
        for genome, coverage in zip(genomes, coverages):
            abundance = coverage / total_coverage
            output.write('\t'.join(map(str, [
                genome.id,
                genome.lineage_string(),
                coverage,
                abundance,
                genome.path
            ])))
            output.write('\n')
            
    @staticmethod
    def generate_scripts(output_folder, simulate_shell_script, genomes, vcovs, sample_id):
        
        single_read_output_folder = "{}/{}/".format(output_folder, "single_reads")

        mkdir_if_not_exists(output_folder)
        mkdir_if_not_exists(single_read_output_folder)

        single_read_output_folder_abs = os.path.abspath(single_read_output_folder)
        joint_read_output_folder_abs = os.path.abspath(output_folder)

        with open(simulate_shell_script, 'w') as output:
            Generator.write_simulate_reads_script(output, joint_read_output_folder_abs, single_read_output_folder_abs, genomes, vcovs, sample_id, gzip=True)
        make_executable(simulate_shell_script)


    @staticmethod
    def generate_meta_file(out, meta_entries, summary=True):

        samples = set(sample.name for _, _, sample in meta_entries)

        if summary:
            out.write("ID\tspecies\tcoverage\tconspecific\tpath\n")
            for sample in samples:
                sp_meta_entries = [me for me in meta_entries if me.sample.name == sample]
                if len(sp_meta_entries) == 0:
                    for genome, vcov, sample in meta_entries:
                        print(genome.to_string(), vcov, sample.name)
                    print(sample)
                    input()
                conspecific_count = len(sp_meta_entries)
                vcov_total = sum(vcov for _, vcov, _ in sp_meta_entries)
                species = sp_meta_entries[0].genome.r_species

                out.write(f"{sample}\t{species}\t{vcov_total}\t{conspecific_count}\t{'None' if conspecific_count > 1 else sp_meta_entries[0].genome.path }\n")

        else:
            out.write("ID\tgenome\tspecies\tcoverage\tpath\n")
            for metaentry in meta_entries:
                genome, vcov, sample = metaentry
                out.write(f"{sample.name}\t{genome.id}\t{genome.r_species}\t{vcov}\t{genome.path}\n")

    @staticmethod
    def generate_meta_files(output_folder, meta_dict):
        for species, metaentries in meta_dict.items():
            meta_file = "{}/{}.meta.tsv".format(output_folder, species.replace(' ', '_'))
            meta_conspecific_file = "{}/{}.conspecific.meta.tsv".format(output_folder, species.replace(' ', '_'))
            with open(meta_file, 'w') as out:
                Generator.generate_meta_file(out, metaentries, summary=True)
            with open(meta_conspecific_file, 'w') as out:
                Generator.generate_meta_file(out, metaentries, summary=False)


    @staticmethod
    def generate_tree_script(out, tree_file, tree_tmp, tree_base, msa):
        out.write("#!/bin/bash\n\n")
        out.write("""
TIMESTR=`date +%Y%m%d-%H%M%S`
TREETMP={}/{}_iqtree_$TIMESTR
MSA={}
iqtree \\
    -s $MSA \\
    -fast \\
    -m {}

mv $MSA.bionj $TREETMP
mv $MSA.ckp.gz $TREETMP
mv $MSA.iqtree $TREETMP
mv $MSA.log $TREETMP
mv $MSA.uniqueseq.phy $TREETMP
cp $MSA.treefile {}
mv $MSA.treefile $TREETMP
        """.format(tree_tmp, tree_base, msa, "GTR", tree_file))

    @staticmethod
    def generate_tree_scripts(script_folder, tree_output_folder, msa_folder, msa_dict, method="iqtree"):
        for species, msa_file in msa_dict.items():
            tree_file = "{}/{}.{}.nwk".format(tree_output_folder, species, method)
            tree_script_file = "{}/{}.{}.sh".format(script_folder, species.replace(' ', '_'), method)
            with open(tree_script_file, 'w') as out:
                Generator.generate_tree_script(out, tree_file,  tree_output_folder, species, "{}/{}/{}".format(msa_folder, species, msa_file))

    @staticmethod
    def generate_roary_scripts(output_folder, script_folder, meta_dict, threads=16, batch_size=10):
        prokka_out = "{}/{}".format(output_folder, "prokka")
        roary_out = "{}/{}".format(output_folder, "roary")

        prokka_script_folder = "{}/{}".format(script_folder, "prokka")
        roary_script_folder = "{}/{}".format(script_folder, "roary")

        mkdir_if_not_exists(prokka_out)
        mkdir_if_not_exists(roary_out)
        mkdir_if_not_exists(prokka_script_folder)
        mkdir_if_not_exists(roary_script_folder)

        all_genomes = list(set(genome for species, meta in meta_dict.items() for genome, vcov, sample in meta))

        for bi, batch in enumerate(batch_list(all_genomes, batch_size)):
            batch_script = "{}/prokka_batch_{}.sh".format(prokka_script_folder, bi)
            with open(batch_script, 'w') as out:
                out.write("#!/bin/bash\n\n")
                out.write("PROKKA_FOLDER={}\n".format(prokka_out))
                out.write("THREADS={}\n".format(threads))
                for index in batch:
                    genome = all_genomes[index]
                    out.write("prokka --outdir $PROKKA_FOLDER --prefix {} --cpus $THREADS --force {}\n".format(
                        genome.id,
                        genome.unzipped_path()
                    ))

        species_to_msa = dict()

        for species, metaentries in meta_dict.items():
            species_name = species.replace(' ', '_')
            species_msa_file = "{}.msa.fna".format(species_name)
            species_script = "{}/{}.sh".format(roary_script_folder, species_name)
            species_roary_out = "{}/{}".format(roary_out, species_name)
            species_genomes_roary_out = "{}/genomes/".format(species_roary_out)
            
            mkdir_if_not_exists(species_roary_out)
            mkdir_if_not_exists(species_genomes_roary_out)

            species_to_msa[species_name] = species_msa_file

            with open(species_script, 'w') as out:
                out.write('#!/bin/bash\n')
                out.write("PROKKA_FOLDER={}\n".format(species_genomes_roary_out))
                out.write("ROARY_FOLDER={}\n".format(species_roary_out))
                for me in metaentries:
                    genome, vcov, sample = me
                    gff_file = "{}.gff".format(genome.id)
                    gff_source = "{}/{}".format(prokka_out, gff_file)
                    out.write("ln -s {} {}\n".format(gff_source, "${{PROKKA_FOLDER}}".format(prokka_out)))

                out.write("\nroary -e --mafft -p 16 ${PROKKA_FOLDER}/*.gff -f ${ROARY_FOLDER}/output\n")
                out.write("ln -s ${{ROARY_FOLDER}}/core_gene_alignment.aln ${{ROARY_FOLDER}}/{}\n".format(species_msa_file))

        return species_to_msa

    @staticmethod
    def generate_benchpro_scripts():
        print("roar")

    @staticmethod
    def resolve_group_vcov(genome_to_group, group_to_vcov, ratio_min=1, ratio_max=9, ratio_type=int):
        genome_to_vcov = dict()
        sample_method = random.randint if ratio_type is int else random.uniform

        for group, group_vcov in set(group_to_vcov.items()):
            genomes = [genome for genome, ggroup in genome_to_group.items() if ggroup == group]

            if len(genomes) == 1:
                genome_to_vcov[genomes[0]] = group_vcov

            logger.info("Group -> vcov: {} -> {}".format(group, group_vcov))
            logger.info("{}: Group size: {}".format(group, len(genomes) - 1))
            ratios = [sample_method(ratio_min, ratio_max) for _ in range(len(genomes) - 1)]
            ratios.append(1)
            random.shuffle(ratios)
            scaler = 1 / sum(ratios)
            logger.info("Ratios: {} -> Sum: {}" .format(ratios, sum(ratios)))
            logger.info("Scaler: {}".format(scaler))
            ratios = [scaler * r for r in ratios]
            logger.info("Scaled ratios: {}" .format(ratios))

            for genome, ratio in zip(genomes, ratios):
                logger.info("Group vcov {}".format(group_vcov))
                genome_to_vcov[genome] = ratio * group_vcov
                logger.info("Scaled vcov {}: {}".format(genome, ratio*group_vcov))

        return genome_to_vcov

            
