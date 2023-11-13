import sys

if len(sys.argv) < 5:
	print("Usage:\npython replace_taxonomy_with.py taxonomy_file genome_col taxonomy_col replacement_taxonomy")
	print("taxonomy_file: tab delimited file")
	print("taxonomy_col: col number of taxonomy string")
	print("genome_col: col number for respective genome")
	print("replacement_taxonomy: tab delimited file genome to taoxnomy")
	exit(9)

taxonomy_path = sys.argv[1]
genome_col = int(sys.argv[2])
taxonomy_col = int(sys.argv[3])
replacement_path = sys.argv[4]

genome2newtax = dict()
with open(replacement_path, 'r') as file:
	for line in file:
		line = line.rstrip()

		genome, taxonomy = line.split('\t')

		genome2newtax[genome] = taxonomy


with open(taxonomy_path, 'r') as file:
	header = True
	for line in file:
		if header:
			header = False
			continue

		line = line.rstrip()

		tokens = line.split('\t')

		genome = tokens[genome_col]
		if genome in genome2newtax:
			new_taxonomy = genome2newtax[genome]
		else:
			print(genome)

		tokens[taxonomy_col] = new_taxonomy

		print('\t'.join(tokens))




