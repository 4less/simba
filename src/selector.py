import os
import random
from loguru import logger
from src.genome_resource import GenomeResource, Rank

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
    def parse_range(range_str: str, condense=True):
        numrange = []
        tokens = range_str.split(',')
        for token in tokens:
            range_token = token.split('-')
            if len(range_token) == 1:
                numrange.append(int(range_token[0]))
            elif len(range_token) == 2:
                numrange.extend(range(int(range_token[0]), int(range_token[1])+1))
            else:
                print("Error")
                exit(9)

        if condense:
            numrange = sorted(list(set(numrange)))

        return numrange

    @staticmethod
    def select_genomes_for_species(genome_resource, species_list, select_conspecific="1", predicate=lambda _: True, predicate_selection=lambda _: True):
        conspecific_pool = Selector.parse_range(select_conspecific)
        
        source_genomes = [g for g in genome_resource.get_genomes() if predicate(g)]

        selected_genomes = set()
        for species in species_list:
            species_genomes = [g for g in source_genomes if g.r_species == species]

            if max(conspecific_pool) >= len(species_genomes):
                logger.warning("Selection range for conspecific strains has options that exceed the number of genomes for the species {} ({})".format(species, len(species_genomes)))

            conspecific_pool_sub = [e for e in conspecific_pool if e <= len(species_genomes)]
            if len(conspecific_pool_sub) == 0:
                logger.error("Error: pool is empty. Skip")
                continue

            select_x = random.choice(conspecific_pool_sub)
            any(selected_genomes.add(g) for g in random.sample(species_genomes, select_x))
            
        return selected_genomes
    
    