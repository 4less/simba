from src.genome_resource import GenomeResource, Rank


class \
SpecificRule:
    def __init__(self):
        self.name = ""

        # select specific taxon and rank
        self.taxon = "s__Escherichia coli"
        self.rank = "species"

        # How many species for this taxon
        self.species_min = None
        self.species_max = None

        # How many conspecific strains per sample
        self.species_conspecific_min = None
        self.species_conspecific_max = None

        # Must be within boundaries of available per species genomes
        self.species_total_genomes_min = None
        self.species_total_genomes_max = None

        self.prevalence_across_samples = None

    def valid(self):
        print(vars(self))
        if any(getattr(self, a) is None for a in vars(self)):
            print('\n'.join(
                ["Attribute {} is not set.".format(a) for a in vars(self) if getattr(self, a) is None]
            ))
            return False
        return True

    def species_rule(self, genome_list):
        species_selection = set(g.r_species for g in genome_list)
        if len(species_selection) < self.species_min:
            print("Not enough species")
            exit(12)
        return True

    def apply(self, genome_resource: GenomeResource, total_samples: int):
        if not self.valid():
            print("not valid")
            exit(12)
        else:
            print("valid")

        taxa = genome_resource.get_taxa_for_rank(Rank.Species)
        eligible_genomes = [g for g in genome_resource.get_genomes() if g.match(self.taxon)]

        # for g in eligible_genomes:
        #     print(g.to_string())


        samples_with = self.prevalence_across_samples * total_samples

        print("Samples with {}: {}/{}".format(self.taxon, samples_with, total_samples))

        self.species_rule

        print("select")





class FillerRule:
    def __init__(self):
        self.name = ""

class RuleSelector:
    def __init__(self):
        self.name = "hello"