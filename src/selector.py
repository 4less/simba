import os
import random
from genome_resource import GenomeResource, Rank

class Selector:
    @staticmethod
    def select_species(genome_resource, x_species=-1, predicate=lambda _: True):
        genomes = genome_resource.get_genomes()
        filtered_genomes = [g for g in genomes if predicate(g)]
        
        if x_species == -1:
            return set(g.r_species for g in filtered_genomes)
        
        species = list(set([g.r_species for g in filtered_genomes]))
        
        selected_species = random.sample(species, x_species)
        
        return selected_species

    @staticmethod
    def select_genomes_for_species(genome_resource, species_list, select_per_species=1, predicate=lambda _: True, predicate_selection=lambda _: True):
        select_x = random.randint(1, select_per_species)
        
        source_genomes = [g for g in genome_resource.get_genomes() if predicate(g)]

        selected_genomes = set()
        for species in species_list:
            species_genomes = [g for g in source_genomes if g.r_species == species]
            any(selected_genomes.add(g) for g in random.sample(species_genomes, select_x))
            
        return selected_genomes
    
    