import os
import numpy as np

from scipy.stats import nbinom
from scipy.stats import powerlaw
from scipy.stats import pareto

import random

from src.simulator_wrappers import Simulator
from src.utils import make_executable, longest_common_prefix, keep_only_folder, mkdir_if_not_exists
from src.genome_resource import Genome


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
        m = len(genomes)
        
        initial_vcov = pareto_distr(m)
        
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
            shortened_genome_path = genome.path.replace(lcpath, "") if len(lcpath) > 0 else genome.path
            genome_str_zip = "${{{}}}/{}".format(VARIABLE_STRING, shortened_genome_path)
            genome_str_unzip = "${{{}}}/{}".format(VARIABLE_STRING, shortened_genome_path).replace('.gz', '')

            shell_command, read_fwd, read_rev = Simulator.get_art_illumina(genome_str_unzip, "${{{}}}/{}".format(SIMULATION_OUTPUT, genome.id), True, 150, coverage)

            read_list_fwd.append(read_fwd)
            read_list_rev.append(read_rev)


            # unzip genome
            output.write("gunzip {}\n".format(genome_str_zip))

            output.write(' '.join(map(str, shell_command)))
            output.write('\n')

            # zip genome again
            output.write("gzip {}\n".format(genome_str_unzip))

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
    def generate_roary_scripts(output_folder, script_folder, species_to_genomes):
        prokka_out = "{}/{}".format(output_folder, "prokka")
        roary_out = "{}/{}".format(output_folder, "roary")

        mkdir_if_not_exists(prokka_out)
        mkdir_if_not_exists(roary_out)

        for species, genomes in species_to_genomes.items():
            species_name = species.replace(' ', '_')
            species_script = "{}/{}.sh".format(script_folder, species_name)
            print(species_script)

            with open(species_script, 'w') as out:
                out.write('#!/bin/bash\n')
                out.write("""
GENOMES=\"{}\"
PROKKA_FOLDER={}
ROARY_FOLDER={}

mkdir -p $PROKKA_FOLDER
mkdir -p $ROARY_FOLDER

for GENOME in $GENOMES; do
	#GENOME=$(echo $line | cut -f2 -d' ')
	GENOME_UNZIP=$(echo $GENOME | sed 's/\.gz//')
    PREFIX=$(basename $GENOME_UNZIP | sed 's/\.fna//')

    echo "gunzip $GENOME"
    echo "srun prokka --outdir $PROKKA_FOLDER --prefix $PREFIX --cpus 16 --force $GENOME_UNZIP"
    echo "gzip $GENOME_UNZIP"

    gunzip $GENOME
	prokka --outdir $PROKKA_FOLDER --prefix $PREFIX --cpus 16 --force $GENOME_UNZIP
    gzip $GENOME_UNZIP
done

roary -e --mafft -p 16 ${{PROKKA_FOLDER}}/*.gff -f ${{ROARY_FOLDER}}
ln -s ${{ROARY_FOLDER}}/*/core_gene_alignment.aln ${{ROARY_FOLDER}}
""".format(' '.join(map(lambda x: x.path, genomes)), prokka_out, roary_out))

    @staticmethod
    def generate_benchpro_scripts():
        print("roar")


            
