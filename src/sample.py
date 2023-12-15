import sys
from src.genome_resource import Rank, Genome
from collections import namedtuple


class Sample:
    def __init__(self):
        self.name = ""
        self.id = 0
        self.filename = ""
        self.genomes: list[Genome] = []
        self.vertical_coverages: list[float] = []
        self.genome2cov = dict()

    def setup(self):
        self.genome2cov = {(genome.id, vcov) for genome, vcov in zip(self.genomes, self.vertical_coverages)}

    def rank_set(self, rank):
        Taxabundance = namedtuple("taxabundance", "genome, vcov")
        
        taxa = set(genome.taxon_at(rank) for genome in self.genomes)
        
        
        return {taxon: [Taxabundance(genome, vcov)
                        for genome, vcov in zip(self.genomes, self.vertical_coverages)
                        if genome.taxon_at(rank) == taxon]
                for taxon in taxa}
