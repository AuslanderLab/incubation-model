#!/usr/bin/env python3
# Imports --------------------------------------------------------------------------------------------------------------
import argparse
import numpy as np
import re
from collections import OrderedDict
import simplejson as json


# Constants ------------------------------------------------------------------------------------------------------------
DESCRIPTION = "Given a genomic alignment and a corresponding tree, calculate the associated mutation rate." \
              "The mutation rate is calculated in a root-agnostic manner, by counting at each branch point the" \
              "transitions between both sides of the branch."


# Classes   ------------------------------------------------------------------------------------------------------------
class GenomeArray:
    """
    Represents a genomic sequence alignment, with multiple options per each location in the genome.
    """
    def __init__(self, sequence=None):
        if sequence is None:
            sequence = ""
        self.__sequence_array = [
            {char} for char in sequence
        ]

    def append(self, pos):
        """Append a position to the array."""
        self.__sequence_array.append(pos)

    def __len__(self):
        return len(self.__sequence_array)

    def __iter__(self):
        for pos in self.__sequence_array:
            yield pos

    def empty(self):
        """Return True if empty, False otherwise."""
        return len(self) == 0

    def merge(self, genome_array):
        """
        Merge two genome arrays to one using frequencies of nucleotides per location.
        :param genome_array: GenomeArray to merge.
        :return: The merged GenomeArray and the number of mutations between the two GenomeArrays.
        """
        assert type(genome_array) == GenomeArray, genome_array

        if genome_array.empty():
            # second array is empty
            # return self, with 0 mutations
            return self, 0
        else:
            # strings should be of equal size to compare mutations
            assert len(self) == len(genome_array)

        mutation_count = 0
        merged_array = GenomeArray()
        for pos_orig, pos_new in zip(self, genome_array):
            if pos_orig == pos_new:
                # the strings are identical at this coordinate, add this to the merged string.
                merged_array.append(pos_new)
            else:
                # the strings are not identical at this coordinate
                # take the intersection if possible
                intersection = pos_orig & pos_new
                if len(intersection) > 0:
                    merged_array.append(intersection)
                else:
                    # there is no intersection, append the array union
                    merged_array.append(pos_orig | pos_new)

                    # don't count gaps as mutations
                    if {'-'} != pos_orig and {'-'} != pos_new:
                        mutation_count += 1

        return merged_array, mutation_count


# Functions ------------------------------------------------------------------------------------------------------------
def proc_fasta(fasta):
    """
    Process aligned Fasta into dict.
    :param fasta: Path to fasta file.
    :return: Dictionary of fasta sequence name -> GenomeArray of sequence.
    """
    ret = OrderedDict()
    seq_name = None
    seq_value = None

    count = 0
    for line in fasta:
        line = line.strip()
        if line.startswith(">"):
            # Reached a new sequence in Fasta, if this is not the first sequence let's save the previous one
            if seq_name is not None:  # This is not the first sequence
                assert seq_value is not None, seq_name

                ret[seq_name] = GenomeArray(seq_value.upper())
            seq_name = line.lstrip(">").split('.')[0]
            count += 1
            seq_value = ""
        else:
            seq_value += line

    # Write last entry
    ret[seq_name] = GenomeArray(seq_value.upper())

    return ret


def proc_newick(newick_str, fdict, sep='**'):
    """
    Convert Newick format string to a list representing tree, where each leaf is represented as an integer, and create a
    dictionary mapping the integer to species names.
    :param newick_str: string of newick tree
    :param fdict: dictionary of aligned sequences
    :param sep: separator (should not appear in the tree string).
    :return: nested list with integers representing the tree, a map from integers to leaf names in the tree
    """
    seq_names = list(fdict.keys())

    # edit sequence names with clear markers
    for seq_name in seq_names:
        newick_str = newick_str.replace(
            seq_name,
            '{' + sep + seq_name + sep + '}'
        )

    # extract all tree names and assert they are the same as the fasta names
    tree_names = [i.split(sep + '}')[0] for i in newick_str.split('{' + sep)[1:]]
    assert set(tree_names) == set(seq_names)

    # create a dictionary matching the index in the tree to the corresponding species
    tree_dict = dict(enumerate(tree_names))

    # remove all extraneous characters from the newick string
    clean_newick_str = re.sub(r'[^(){}*,]', '', newick_str)

    # replace loci where species names used to occur with a count
    count = 0
    while clean_newick_str.find('{'+sep+sep+'}') > -1:
        clean_newick_str = clean_newick_str.replace('{'+sep+sep+'}', str(count), 1)
        count += 1

    # convert the string representation of a nested list of integers to a python list
    newick_list = json.loads(clean_newick_str.replace('(', '[').replace(')', ']'))

    return newick_list, tree_dict


def extract_seq(obj, tree_dict, fdict):
    """Extract the GenomeArray from the corresponding integer. If already a GenomeArray, return it."""
    if isinstance(obj, int):
        return fdict[tree_dict[obj]]
    else:
        assert isinstance(obj, GenomeArray)
        return obj


def is_leaf(subtree):
    """Return True if a subtree is a leaf node."""
    return isinstance(subtree, int)


def is_sequence(subtree):
    """Return True if a subtree is a sequence"""
    return isinstance(subtree, GenomeArray)


def is_subtree(subtree):
    """Return True if a subtree is a subtree."""
    return isinstance(subtree, list)


def is_term_node(subtree):
    """Return True if subtree is a terminal node, either a leaf node or the corresponding sequence."""
    return is_sequence(subtree) or is_leaf(subtree)


def calc_subtree_mutation_rate(tree_list, mutation_counts, tree_dict, fdict):
    """
    Recursively extract mutation rates from nested list representation of tree.
    :param tree_list: Nested list of integers that represent sequences in the tree.
    :param mutation_counts: list of mutation rates
    :param tree_dict: dictionary mapping integers to names
    :param fdict: dictionary mapping names to aligned strArray sequence
    :return: average mutation rate, list of mutation rates and consensus ancestor string
    """
    first_subtree, second_subtree = tree_list[0],  tree_list[1]

    if is_term_node(first_subtree) and is_term_node(second_subtree):
        # both subtrees are terminal nodes - extract sequences if needed, and merge
        first_sequence = extract_seq(first_subtree, tree_dict, fdict)
        second_sequence = extract_seq(second_subtree, tree_dict, fdict)
        merged_array, mutation_count = first_sequence.merge(second_sequence)
        mutation_counts.append(mutation_count)
        ret = mutation_counts, merged_array

    elif ((is_term_node(first_subtree) and is_subtree(second_subtree)) or
            (is_term_node(second_subtree) and is_subtree(first_subtree))):
        # one subtree is a terminal node, the other is a subtree that continues to branch
        subtree, node = (first_subtree, second_subtree) \
            if is_subtree(first_subtree) \
            else (second_subtree, first_subtree)

        # get mutation counts and merged string for the subtree
        mutation_counts, merged_array = calc_subtree_mutation_rate(subtree, mutation_counts, tree_dict, fdict)

        # get the sequence from the terminal node
        sequence = extract_seq(node, tree_dict, fdict)

        # merge the sequences and add the new mutation count
        merged_array, mutation_count = sequence.merge(merged_array)
        mutation_counts.append(mutation_count)

        # return the values
        ret = mutation_counts, merged_array

    elif (is_subtree(first_subtree)) and (is_subtree(second_subtree)):
        # both are subtrees
        first_mutation_count, first_merged_string = calc_subtree_mutation_rate(
            first_subtree, mutation_counts, tree_dict, fdict
        )
        second_mutation_count, second_merged_string = calc_subtree_mutation_rate(
            second_subtree, first_mutation_count, tree_dict, fdict
        )

        # merge the subtrees
        merged_array, mutation_count = first_merged_string.merge(second_merged_string)
        second_mutation_count.append(mutation_count)

        # return the values
        ret = second_mutation_count, merged_array

    return ret


def calc_mutation_rate(fasta_path, newick_path):
    """
    Calculates mutation rates from fasta alignment and respective newick tree.
    :param fasta_path: fasta alignment path and respective newick tree path
    :param newick_path: corresponding newick tree
    :return: average mutation rate, list of mutation rates and consensus ancestor string.
    """
    # read fasta into dict
    with open(fasta_path) as handle:
        fdict = proc_fasta(handle)

    # read newick tree
    with open(newick_path) as handle:
        tree_lines = handle.readlines()
        assert len(tree_lines) == 1, tree_lines
        tree_text = tree_lines[0]

    newick_list, tree_dict = proc_newick(tree_text, fdict)

    # traverse the tree
    merged_array = GenomeArray()
    mutation_counts = []
    for sub_tree in newick_list:
        if is_leaf(sub_tree):
            # subtree is a leaf
            # recode as two leaves to match sub tree format
            sub_tree = [sub_tree, sub_tree]

        # get mutation rate in subtree
        sub_mutation_counts, sub_merged_array = calc_subtree_mutation_rate(sub_tree, [], tree_dict, fdict)
        mutation_counts.extend(sub_mutation_counts)

        # merge with overall mutation rate, merged_string
        merged_array, new_mutation_counts = sub_merged_array.merge(merged_array)
        mutation_counts.append(new_mutation_counts)

    return np.mean(mutation_counts) / len(merged_array), mutation_counts, merged_array


# Main -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--alignment", type=str, help="Alignment of genomes in fasta format", required=True
    )
    parser.add_argument(
        "--tree", type=str, help="Tree corresponding to genome alignment in newick format", required=True
    )
    args = parser.parse_args()

    mutation_rate = calc_mutation_rate(args.alignment, args.tree)[0]
    print(mutation_rate)
