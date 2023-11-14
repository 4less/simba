import glob
from pathlib import Path

class Genome:
    def __init__(self, genome, path, size):
        self.genome = genome
        self.path = path
        self.size = size
        
    def to_string(self):
        return f"{self.genome}\t{self.path}\t{self.size}"
        
def load_genomes(genome_folder, genomes_file):
    # Assume that in this simple setting the filename without the extension is the name
    # and also the identifier in the "genome" column in the taxonomy file.
    # If there is no taxonomy file then this does not matter
    fna_globber = f"{genome_folder}/*.fna"
    fasta_globber = f"{genome_folder}/*.fasta"
    fa_globber = f"{genome_folder}/*.fa"
    fnas = glob.glob(fna_globber)
    fasta = glob.glob(fasta_globber)
    fa = glob.glob(fa_globber)
    genomes = [Genome(Path(g).stem, g, 0) for g in fnas + fasta + fa]
    return genomes
