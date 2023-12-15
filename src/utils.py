from src.tree_node import TreeNode
import random
import numpy as np
import os
import matplotlib.pyplot as plt

from scipy.stats import nbinom
from scipy.stats import powerlaw
from scipy.stats import pareto

def load_dict(path: str, key_col=0, value_col=1, delimiter='\t', transform_key=lambda x: x, transform_value=lambda x: x):
    try:
        output_dict = dict()
        with open(path, 'r') as file:
            for line in file:
                line = line.rstrip()

                tokens = line.split(delimiter)
    
                output_dict[transform_key(tokens[key_col])] = transform_value(tokens[key_col])
                
        return output_dict
    except:
        print("Unable to to open file {}".format(path))


def Sample(set_var):
    random.choice(tuple(set_var))


def BuildNewickString(node: TreeNode):
    if len(node.children) == 0:
        return "\"{}\"".format(node.name)
    else:
        return "({}){}".format(
            ','.join([BuildNewickString(child) for child in node.children]),
            "\"{}\"".format(node.name))


def NewickTree(genome_list, taxonomy, tree_nodes):
    global lines
    lines = []

    root_node = TreeNode("r__Root")
    tree_nodes["r__Root"] = root_node

    lineage = ""
    node = None

    for genome in genome_list:
        lineage = taxonomy.PrintLineage(genome)

        tokens = lineage.split(';')

        node = None
        parent_node = root_node
        for taxon in tokens:
            if taxon not in tree_nodes:
                tree_nodes[taxon] = TreeNode(taxon)

            node = tree_nodes[taxon]

            parent_node.AddChild(node)
            parent_node = node

    return "({});".format(BuildNewickString(root_node))


def plot_distr(x, y, output: str):
    plt.plot(x, y, label='randomized')
    plt.xlabel('samples', fontsize=15)
    plt.ylabel('PDF', fontsize=15)
    plt.title('Distribution', fontsize=15)
    plt.grid(b=True, color='grey', alpha=0.3, linestyle='-.', linewidth=2)
    plt.rcParams["figure.figsize"] = [5, 5]
    plt.legend(loc='best')
    plt.savefig(output)


def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2  # copy R bits to X
    os.chmod(path, mode)


def fasta_length(path: str):
    length = 0
    with open(path, 'r') as file:
        for line in file:
            line = line.rstrip()

            if not line.startswith('>'):
                length += len(line)
    return length


def read_count(rel_abundance: float, read_length: int, genome_length: int):
    return (genome_length * rel_abundance) / read_length


def vertical_coverage(read_count: int, read_length: int, genome_length: int):
    return read_count * read_length / genome_length


def upscale_read_count(counts, new_total: int):
    factor: float = new_total / sum(counts)
    return [factor * count for count in counts]


def mkdir_if_not_exists(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)


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
    output = np.array([pareto.pdf(x=samples, b=a, loc=0, scale=x_m)])

    values = randomize_distr(output.T[1:])
    scale = 1 / sum(values)
    values = [scale * val for val in values]

    return range(0, data_points), values


def longest_common_prefix(strings):
    if not strings:
        return ""

    # Sort the strings to ensure the shortest string is first
    strings.sort()

    # Consider the first and last strings after sorting
    first_str = strings[0]
    last_str = strings[-1]

    # Find the common prefix between the first and last strings
    prefix = ""
    for i in range(len(first_str)):
        if first_str[i] == last_str[i]:
            prefix += first_str[i]
        else:
            break

    return prefix


def keep_only_folder(path):
    index_of_last_slash = path.rfind('/')
    modified_path = path
    if index_of_last_slash != -1:
        modified_path = path[:index_of_last_slash]
    return modified_path
