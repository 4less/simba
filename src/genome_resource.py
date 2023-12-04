import os
import sys
import pandas as pd
from collections import Counter


class Genome:
    RANK_FIELDS = ['r_domain', 'r_phylum', 'r_class', 'r_order', 'r_family', 'r_genus', 'r_species']
    
    def __init__(self, gid, path, r_domain, r_phylum, r_class, r_order, r_family, r_genus, r_species):
        self.id = gid
        self.path = path
        self.r_domain = r_domain
        self.r_phylum = r_phylum
        self.r_class = r_class
        self.r_order = r_order
        self.r_family = r_family
        self.r_genus = r_genus
        self.r_species = r_species


    @staticmethod
    def from_lineage_string(gid, path, lineage):
        tokens = lineage.split(';')
        if len(tokens) != 7:
            print(f"lineage string {lineage} is malformated")
        return Genome(gid, path, tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5], tokens[6])

    def lineage_string(self):
        return ';'.join([getattr(self, f) for f in self.RANK_FIELDS])

    def match(self, taxon):
        return any(t == taxon for t in [getattr(self, f) for f in self.RANK_FIELDS])

    def to_string(self):
        return "{}\t{}".format(self.id, self.lineage_string())


class Rank:
    Domain = "Domain"
    Phylum = "Phylum"
    Class = "Class"
    Order = "Order"
    Family = "Family"
    Genus = "Genus"
    Species = "Species"
    All = [Domain, Phylum, Class, Order, Family, Genus, Species]


class GenomeResource:
    GENOME_COLUMN_NAME = ""
    GENOMEID = 'ID'
    GENOMEPATH = 'Path'
    GENOME_TAXONOMY = 'Taxonomy'
    RANKS = Rank.All

    def __init__(self, path: str):
        self.map_delim = '\t'
        self.tax_delim = ';'

        self.df = None
        self.load(path)

    def load(self, path: str):
        # Load map file
        self.df = pd.read_csv(path, sep=self.map_delim)

        self.df[self.RANKS] = self.df[self.GENOME_TAXONOMY].str.split(self.tax_delim, expand=True)

    def print_summary(self):
        for column_name, count in self.df[self.RANKS].nunique().items():
            print(f"{column_name} has {count} unique values.")
        for rank in self.RANKS:
            print("Entries for {}: {}".format(rank, self.df[rank].value_counts()))

    def x_of_rank(self, rank, exceptions):
        return len(self.df[rank].tolist())

    def has_x_of_rank(self, x, rank, exceptions=None):
        exception_set = set() if exceptions is None else exceptions
        return len(set(self.df[rank].tolist())) >= x

    def get_taxa_for_rank(self, rank):
        return set(self.df[rank].tolist())

    def get_genome_ids(self):
        return self.df[self.GENOMEID].tolist()

    def get_genomes(self):
        return [
            Genome.from_lineage_string(
                row[self.GENOMEID],
                row[self.GENOMEPATH],
                row[self.GENOME_TAXONOMY]) for index, row in self.df.iterrows()]
