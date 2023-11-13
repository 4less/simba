import sys
from ete3 import NCBITaxa

genome_col = 0
abundance_col = 2
lineage_col = 5

ncbiTaxa = NCBITaxa()

class Ranks:
    ncbi = ["superkingdom","phylum","class","order","family","genus","species", "strain"]
    gtdb = ["d__", "p__", "c__", "o__", "f__", "g__", "s__", "genome"]
    ncbi2gtdb = {
        'root': 'root',
        'superkingdom': 'd__',
        'kingdom': 'k__',
        'phylum': 'p__',
        'class': 'c__',
        'order': 'o__',
        'family': 'f__',
        'genus': 'g__',
        'species': 's__',
        'strain': 'genome'
    }
    gtdb2ncbi = {
        'root': 'root',
        'd__': 'superkingdom',
        'k__': 'kingdom',
        'p__': 'phylum',
        'c__': 'class',
        'o__': 'order',
        'f__': 'family',
        'g__': 'genus',
        's__': 'species',
        'genome': 'strain'
    }

    def __init__(self):
        pass


#genome	size	rel_abundance	read_count	sequencing_depth	lineage	absolute_path


composition_file = sys.argv[1]

rank_dict = dict()

with open(composition_file, 'r') as file:
	for rank in Ranks.gtdb[:-1]:
		for line in file:
			line = line.rstrip()

			tokens = line.split('\t')

			genome = tokens[genome_col]
			abundance = float(tokens[abundance_col])
			taxonomy_str = tokens[lineage_col]
			taxonomy = [taxon for taxon in taxonomy_str.split(';') if len(taxon) > 3]

			for x in range(len(taxonomy)):
				if taxonomy[x].startswith(rank):
					
					break

					
		file.seek(0)