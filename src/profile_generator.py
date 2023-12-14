import os
import numpy as np
import re

from scipy.stats import nbinom
from scipy.stats import powerlaw
from scipy.stats import pareto

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

def randomize_distr(distr, alpha: int = 2):
    sum_old = 0
    sum_new = 0
    for i in range(0, len(distr)):
        # weight = 1 / (i/2+1)
        # weight = 3/(i+1)
        weight = 0.3 + (10 / len(distr) * min(i, 10)) / 10
        random_element = random.uniform(-1 * weight, weight)
        distr[i] += abs(random_element * distr[i] * alpha)
        distr = list(distr)
        distr.sort(reverse=True)

    return distr

def linear_distr(data_points: int):
    return list(range(1, data_points + 1))

def pareto_distr(data_points: int):
    a = 0.15
    x_m = 1
    samples = np.linspace(start=0, stop=data_points + 1, num=data_points + 1)
    output = np.array(pareto.pdf(x=samples, b=a, loc=0, scale=x_m))
    
    values = randomize_distr(output.T[1:])
    scale = 1 / sum(values)
    
    values = [scale * val for val in values]

    return values


# def scale_counts(counts, new_total: int):
#     factor: float = new_total / sum(counts)
#     return [factor * count for count in counts]

def scale_counts(counts, min_vcov, max_vcov):
    min_c = min(counts)
    counts = [c - min_c for c in counts]
    max_c = max(counts)
    
    factor: float = (max_vcov - min_vcov) / max_c
    
    return [(factor * count) + min_vcov for count in counts]


def batch_list(input_list, batch_size):
    return [list(range(i, min(i + batch_size, len(input_list)))) for i in range(0, len(input_list), batch_size)]

class ProfileGenerator:
    @staticmethod
    def generate_vertical_coverage(genomes, min_vcov, max_vcov, method="pareto"):
        """Generates vertical coverages for genomes

        Parameters
        ----------
        genomes : str
            The file location of the spreadsheet
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
        list
            a list of strings used that are the header columns
        """

        genomes = list(genomes)
        random.shuffle(genomes)
        m = len(genomes)
        
        if method == "pareto":
            initial_vcov = pareto_distr(m)
        elif method == "linear":
            initial_vcov = linear_distr(m)
        else:
            initial_vcov = linear_distr(m)
        
        scaled_vcov = scale_counts(initial_vcov, min_vcov, max_vcov)
        
        return genomes, scaled_vcov
        

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
            ProfileGenerator.write_simulate_reads_script(output, joint_read_output_folder_abs, single_read_output_folder_abs, genomes, vcovs, sample_id, gzip=True)
        make_executable(simulate_shell_script)

    @staticmethod
    def generate_roary_scripts(output_folder, script_folder, species_to_genomes, threads=16, batch_size=10):
        prokka_out = "{}/{}".format(output_folder, "prokka")
        roary_out = "{}/{}".format(output_folder, "roary")

        prokka_script_folder = "{}/{}".format(script_folder, "prokka")
        roary_script_folder = "{}/{}".format(script_folder, "roary")

        mkdir_if_not_exists(prokka_out)
        mkdir_if_not_exists(roary_out)
        mkdir_if_not_exists(prokka_script_folder)
        mkdir_if_not_exists(roary_script_folder)

        all_genomes = list(set(genome for species, genomes in species_to_genomes.items() for genome in genomes))

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

        for species, genomes in species_to_genomes.items():
            species_name = species.replace(' ', '_')
            species_script = "{}/{}.sh".format(roary_script_folder, species_name)
            species_roary_out = "{}/{}".format(roary_out, species_name)
            species_genomes_roary_out = "{}/genomes/".format(species_roary_out)
            
            mkdir_if_not_exists(species_roary_out)
            mkdir_if_not_exists(species_genomes_roary_out)

            with open(species_script, 'w') as out:
                out.write('#!/bin/bash\n')
                out.write("PROKKA_FOLDER={}\n".format(species_genomes_roary_out))
                out.write("ROARY_FOLDER={}\n".format(species_roary_out))
                for genome in genomes:
                    gff_file = "{}.gff".format(genome.id)
                    gff_source = "{}/{}".format(prokka_out, gff_file)
                    out.write("ln -s {} {}\n".format(gff_source, "${{PROKKA_FOLDER}}".format(prokka_out)))

                out.write("\nroary -e --mafft -p 16 ${PROKKA_FOLDER}/*.gff -f ${ROARY_FOLDER}")


    @staticmethod
    def generate_benchpro_scripts():
        print("roar")


            
