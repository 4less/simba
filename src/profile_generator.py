import os
import numpy as np

from scipy.stats import nbinom
from scipy.stats import powerlaw
from scipy.stats import pareto

import random

from src.simulator_wrappers import Simulator
from src.utils import make_executable, longest_common_prefix, keep_only_folder, mkdir_if_not_exists


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
    def generate_vertical_coverage(genomes, min_vcov, max_vcov):
        genomes = list(genomes)
        m = len(genomes)
        
        initial_vcov = pareto_distr(m)
        
        scaled_vcov = scale_counts(initial_vcov, min_vcov, max_vcov)
        
        return genomes, scaled_vcov
        

    @staticmethod
    def write_simulate_reads_script(output, simulation_output, genomes, coverages, sample_id, gzip=False):
        output.write('#!/bin/bash\n\n')
        
        genomes = list(genomes)
        
        paths = [g.path for g in genomes]
        lcp = longest_common_prefix(paths)
        lcpath = keep_only_folder(lcp)
        lcpath = lcpath + '/' if len(lcpath) > 0 else ""
        
        VARIABLE_STRING="PATH_PREFIX"
        output.write("{}={}\n".format(VARIABLE_STRING, lcpath))
        
        
        SIMULATION_OUTPUT="OUTPUT_FOLDER"
        output.write("{}={}\n".format(SIMULATION_OUTPUT, simulation_output))
        
        read_list_fwd = []
        read_list_rev = []

        for genome, coverage in zip(genomes, coverages):
            shortened_path = genome.path.replace(lcpath, "") if len(lcpath) > 0 else genome.path
            shell_command, read_fwd, read_rev = Simulator.get_art_illumina("${{{}}}/{}".format(VARIABLE_STRING, shortened_path), "${{{}}}/{}".format(SIMULATION_OUTPUT, genome.id), True, 150, coverage)

            read_list_fwd.append(read_fwd)
            read_list_rev.append(read_rev)

            output.write(' '.join(map(str, shell_command)))
            output.write('\n')

        read_fwd_out = "{}_1.fq".format(sample_id)
        read_rev_out = "{}_2.fq".format(sample_id)

        output.write("cat {} > {}\n".format(' '.join(read_list_fwd), read_fwd_out))
        output.write("cat {} > {}\n".format(' '.join(read_list_rev), read_rev_out))

        if gzip:
            output.write("gzip {}\n".format(' '.join(read_list_fwd)))
            output.write("gzip {}\n".format(' '.join(read_list_rev)))
            output.write("gzip {}\n".format(read_fwd_out))
            output.write("gzip {}\n".format(read_rev_out))

        output.write('\n')




    @staticmethod
    def write_profile(output, genomes, coverages, delimiter='\t'):
        for genome, coverage in zip(genomes, coverages):
            output.write(f"{genome.id}{delimiter}{genome.lineage_string()}{delimiter}{coverage}\n")
            
    @staticmethod
    def generate(output_folder, simulate_shell_script, genomes, vcovs, sample_id):
        gold_standard_profile = output_folder + '/gold_standard_profile.tsv'
        
        read_output_path = "{}/{}/".format(output_folder, "reads")
        mkdir_if_not_exists(read_output_path)
        
        with open(simulate_shell_script, 'w') as output:
            ProfileGenerator.write_simulate_reads_script(output, read_output_path, genomes, vcovs, sample_id, gzip=True)
        make_executable(simulate_shell_script)

        with open(gold_standard_profile, 'w') as output:
            ProfileGenerator.write_profile(output, genomes, vcovs)
            
