import argparse
import os

def get_simba_argument_parser():
    my_parser = argparse.ArgumentParser(description='Select random genomes')
    
    ##############################################################################
    # I/O
    ##############################################################################
    
    # Input folder
    my_parser.add_argument(
        '-f', 
        '--genome_folder', 
        type=str,
        required=True,
        action='store',
        help='Genome folder')
    
    # Genome list
    my_parser.add_argument(
        '-l', 
        '--genome_list', 
        type=str,
        required=False, # changed from True
        action='store',
        help='Genome list (List of genomes to pick from)')
    
    # Genome list
    my_parser.add_argument(
        '-y', 
        '--selected_genomes', 
        type=str,
        required=False,
        action='store',
        help='Preselected genomes')
    
    # Instructions
    my_parser.add_argument(
        '-i', 
        '--instructions', 
        type=str,
        required=False,
        action='store',
        help='Instructions for how to simulate genomes')
    
    # Instructions
    my_parser.add_argument(
        '-m', 
        '--samples', 
        type=int,
        required=False,
        action='store',
        help='Number of samples to simulate.')
    
    # Taxonomy
    my_parser.add_argument(
        '-t', 
        '--taxonomy', 
        required=True,
        type=str,
        action='store',
        help='Tab separated file with 2 columns. Map genome (col: 1) to GTDB string (col: 2).')
    
    # Output folder
    my_parser.add_argument(
        '-o', 
        '--output', 
        type=str,
        required=True,
        action='store',
        help='Output folder (will be created)')
    
    ##############################################################################
    # Taxonomy params
    ##############################################################################
    
    my_parser.add_argument(
        '-n', 
        '--num_genomes', 
        type=int,
        required=False,
        action='store',
        help='Number of genomes to pick',
        default=10)
    
    my_parser.add_argument(
        '-s', 
        '--max_genomes_per_species', 
        required=False,
        type=int,
        action='store',
        help='Maximum number of genomes per species',
        default=1)
    
    my_parser.add_argument(
        '-g', 
        '--max_genomes_per_genus', 
        required=False,
        type=int,
        action='store',
        help='Maximum number of genomes per genus',
        default=4)
    
    my_parser.add_argument(
        '-q', 
        '--max_genomes_per_family', 
        required=False,
        type=int,
        action='store',
        help='Maximum number of genomes per family',
        default=5)
    
    my_parser.add_argument(
        '-p', 
        '--min_phylum_coverage', 
        required=False,
        type=int,
        action='store',
        help='Minimum number of phyla',
        default=6)
    
    ##############################################################################
    # Simulation params
    ##############################################################################
    
    my_parser.add_argument(
        '-r', 
        '--read_length',
        type=int,
        action='store',
        help='Read length of simulated genomes',
        default=150)
    
    my_parser.add_argument(
        '-v', 
        '--target_coverage',
        type=float,
        action='store',
        help='Target mean horizontal coverage.',
        default=5)
    
    my_parser.add_argument(
        '-e', 
        '--error_rate',
        type=float,
        action='store',
        help='Error rate for simulation',
        default=0.001)
    
    
    my_parser.add_argument(
        '-R', 
        '--indel_rate',
        type=float,
        action='store',
        help='Indel rate for simulation',
        default=0.0001)
    
    my_parser.add_argument(
        '-M', 
        '--mutation_rate',
        type=float,
        action='store',
        help='Mutation rate for simulation',
        default=0)
    
    
    my_parser.add_argument(
        '-x', 
        '--execute',
        action='store_true',
        help='Simulate reads after script.')
    
    my_parser.add_argument(
        '-P', 
        '--no_plots',
        action='store_true',
        help='No plots (commandline)')
    
    my_parser.add_argument(
        '-c', 
        '--copy_genomes',
        action='store_true',
        help='Copy genomes into genomes folder')
    
    return my_parser

# ParserArguments
class Args:
    GENOME_MAP = 'genomes_map'
    CONSPECIFIC = 'conspecific'
    SAMPLES = 'samples'
    SPECIES_NUMBER = 'species_number'
    SHARED_SPECIES_NUMBER_TOTAL = 'shared_species_number_total'
    SHARED_SPECIES_NUMBER_PAIRWISE = 'shared_species_number_pairwise'
    TREE = 'tree'
    OUTPUT_FOLDER = 'output_folder'


def get_simba_refactor_argument_parser():
    my_parser = argparse.ArgumentParser(description='Select random genomes')

    ##############################################################################
    # I/O
    ##############################################################################

    # Genome list
    my_parser.add_argument(
        '-g',
        f"--{Args.GENOME_MAP}",
        type=str,
        required=False, # changed from True
        action='store',
        help='tab-delimited file containing genome id, genome path and taxonomy for each genome. (--genomes_map_help for more information)')

    my_parser.add_argument(
        f"--{Args.TREE}",
        type=str,
        required=False, # changed from True
        action='store',
        help='Supply tree to get more information per sample.')

    my_parser.add_argument(
        f"--{Args.OUTPUT_FOLDER}",
        type=str,
        required=False, # changed from True
        action='store',
        help='Supply tree to get more information per sample.')

    ##############################################################################
    # Within sample parameters
    ##############################################################################

    my_parser.add_argument(f"--{Args.CONSPECIFIC}",
                           action='store_true',
                           help='Enable having conspecific strains in a sample')


    my_parser.add_argument(f"--{Args.SPECIES_NUMBER}",
                           type=int,
                           default=1,
                           action='store',
                           help='Number of different species in a sample')

    ##############################################################################
    # Between sample parameters
    ##############################################################################

    my_parser.add_argument(f"--{Args.SAMPLES}",
                           type=int,
                           default=1,
                           action='store',
                           help='Number of samples')

    my_parser.add_argument(f"--{Args.SHARED_SPECIES_NUMBER_TOTAL}",
                           type=int,
                           default=None,
                           action='store',
                           help='Number of different species in a sample')

    my_parser.add_argument(f"--{Args.SHARED_SPECIES_NUMBER_PAIRWISE}",
                           type=int,
                           default=None,
                           action='store',
                           help='Number of different species in a sample')


    return my_parser

