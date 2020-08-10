#!/usr/bin/env python3
# Imports --------------------------------------------------------------------------------------------------------------
import numpy as np
import random
from Bio.SeqUtils import GC
import argparse
from Bio import SeqIO
import warnings
from collections import OrderedDict
from Bio.SeqUtils.CodonUsage import SynonymousCodons, CodonAdaptationIndex

# Constants ------------------------------------------------------------------------------------------------------------
DESCRIPTION = "Calculate the predicted upper-limit incubation time for a given ssRNA virus." \
              "The Genbank and alignment should all have the an entry with the genome ID."

MODEL_PARAMS = {
    'Coronaviridae': {
        "coef": [1.88026468, 1.17517605, -0.82794489, 0.06439647],
        "intercept": 7.571428571428568,
        "scaler_mean": [9.00000000e+00, 3.73629557e+01, 3.45446014e-05, 6.56935391e-01],
        "scaler_var": [5.14285714e+00, 9.23125265e+00, 1.72812845e-09, 2.22493344e-04]
    },

    'Coronaviridae_Paramyxoviridae': {
        "coef": [1.06089815,  1.35295837, -0.87588244,  0.63610275],
        "intercept": 7.727272727272722,
        "scaler_mean": [8.45454545e+00, 3.81036647e+01, 3.66221062e-05, 6.62414995e-01],
        "scaler_var": [4.79338843e+00, 1.53661941e+01, 1.56856212e-09, 3.20114171e-04]
    },

    'Coronaviridae_Pneumoviridae':{
        "coef": [1.56944687,  1.14647992, -0.7224304, 0.13611629],
        "intercept": 7.444444444444447,
        "scaler_mean": [9.22222222e+00, 3.67818551e+01, 2.98131326e-05, 6.57218909e-01],
        "scaler_var": [5.06172840e+00, 8.89159026e+00, 1.45909879e-09, 1.75575650e-04]
    },


    'Coronaviridae_Low':{
        "coef": [0.227763991338773,  0.061712377887596495, -0.0, 0.],
        "intercept": 3.857142857142857,
        "scaler_mean": [9.00000000e+00, 3.73629557e+01, 3.45446014e-05, 6.56935391e-01],
        "scaler_var": [5.14285714e+00, 9.23125265e+00, 1.72812845e-09, 2.22493344e-04]
    },

    'Coronaviridae_High':{
        "coef": [1.6258975708554884,  0.945883434848541, -0.6393060388885797, 0.07670065125482098],
        "intercept": 8.857142857142854,
        "scaler_mean": [9.00000000e+00, 3.73629557e+01, 3.45446014e-05, 6.56935391e-01],
        "scaler_var": [5.14285714e+00, 9.23125265e+00, 1.72812845e-09, 2.22493344e-04]
    }
}

COLUMNS = [
    "num_genes",
    "gc_content",
    "Position_change_var",
    "cai",
]

NUCLEOTIDES = {"A", "C", "G", "T", "U", "-"}

CAI_FREQS = {
    "TTT": 17.6, "TCT": 15.2, "TAT": 12.2, "TGT": 10.6, "TTC": 20.3, "TCC": 17.7, "TAC": 15.3,
    "TGC": 12.6, "TTA": 7.7, "TCA": 12.2, "TAA": 1.0, "TGA": 1.6, "TTG": 12.9, "TCG": 4.4,
    "TAG": 0.8, "TGG": 13.2, "CTT": 13.2, "CCT": 17.5, "CAT": 10.9, "CGT": 4.5, "CTC": 19.6,
    "CCC": 19.8, "CAC": 15.1, "CGC": 10.4, "CTA": 7.2, "CCA": 16.9, "CAA": 12.3, "CGA": 6.2,
    "CTG": 39.6, "CCG": 6.9, "CAG": 34.2, "CGG": 11.4, "ATT": 16.0, "ACT": 13.1, "AAT": 17.0,
    "AGT": 12.1, "ATC": 20.8, "ACC": 18.9, "AAC": 19.1, "AGC": 19.5, "ATA": 7.5, "ACA": 15.1,
    "AAA": 24.4, "AGA": 12.2, "ATG": 22.0, "ACG": 6.1, "AAG": 31.9, "AGG": 12.0, "GTT": 11.0,
    "GCT": 18.4, "GAT": 21.8, "GGT": 10.8, "GTC": 14.5, "GCC": 27.7, "GAC": 25.1, "GGC": 22.2, "GTA": 7.1,
    "GCA": 15.8, "GAA": 29.0, "GGA": 16.5, "GTG": 28.1, "GCG": 7.4, "GAG": 39.6, "GGG": 16.5
}


# Functions ------------------------------------------------------------------------------------------------------------
def get_gene_count(genbank):
    """Return gene count from Genbank file"""
    return len([x for x in genbank.features if x.type == "gene"])


def read_gb(path, genome_id):
    """Read Genbank into dictionary of records"""
    gb_records = {x.id: x for x in SeqIO.parse(path, "genbank")}
    assert genome_id in list(gb_records.keys()), "{} must be in genbank file".format(genome_id)
    return gb_records[genome_id]


def predict(model_data, model_params):
    """Predict with a regression model."""
    pred = sum(model_data * model_params["coef"]) + model_params["intercept"]
    return pred


def random_NN(nns):
    """
    Generate a random base.
    :return: Random base.
    """
    return random.choice(nns)


def proc_sequence(sequence):
    """
    Taken from Seeker and converted
    """
    aa_dict = {
        'W': ['A', 'T'],
        'Y': ['C', 'T']
    }

    ret = sequence
    for item in aa_dict:
        ret = ret.replace(item, random_NN(aa_dict[item]))

    return ret.upper()


def is_nucleotide(sequence):
    """Return true if sequence is nucleotide"""
    count = 0
    for nucl in NUCLEOTIDES:
        count += sequence.count(nucl)
    return (count / len(sequence)) > 0.90


def proc_fasta(fasta, alphabet="amino_acid"):
    """Process aligned Fasta into dict."""
    assert alphabet in {"amino_acid", "nucleotide"}, "alphabet must be 'amino' or 'nucl'"

    ret = OrderedDict()
    for entry in SeqIO.parse(fasta, "fasta"):
        entry_sequence = proc_sequence(str(entry.seq))
        seq_is_nuc = is_nucleotide(entry_sequence)
        if alphabet == "amino_acid":
            assert not seq_is_nuc, "File {} should be {}".format(fasta.name, alphabet)
        else:
            assert seq_is_nuc, "File {} should be {}".format(fasta.name, alphabet)

        ret[entry.id] = entry_sequence

    assert len(ret) > 0, "{} is is empty".format(fasta.name)
    return ret


def get_position_variance(aln_seqdict):
    """Calculates estimated variance of mutation counts per family."""
    mafs = []
    sequences = list(aln_seqdict.values())
    for base_idx in range(len(sequences[0])):
        bases_col = [sequence[base_idx] for sequence in sequences]

        base_counts = [
            sum([base == unique_base for base in bases_col]) for unique_base in np.unique(bases_col).tolist()
        ]

        base_counts.pop(base_counts.index(max(base_counts)))
        mafs.append(len(np.unique(base_counts)) / len(sequences))

    return np.var(mafs)


def calc_cai(sequence, genbank, cai_freqs=CAI_FREQS):
    """Return the CAI for a given genome."""
    # create CAI index
    cai_index = {}
    for codons in SynonymousCodons.values():
        codons = list(codons)
        codon_freqs = np.array([cai_freqs[x] for x in codons])
        max_freq = max(codon_freqs)
        codon_freqs = codon_freqs / max_freq
        for i, x in enumerate(codons):
            cai_index[x] = codon_freqs[i]
    cai_table = CodonAdaptationIndex()
    cai_table.set_cai_index(cai_index)

    # concatenate ORFs
    orfs = [x for x in genbank.features if x.type.lower() == "cds"]
    cds_seq = ""
    for orf in orfs:
        cds_seq += proc_sequence(str(orf.extract(sequence)).upper().replace("U", "T"))

    # return cai
    return cai_table.cai_for_gene(cds_seq)


def create_model_data(genome_id, aln_seqdict, genbank, model_params):
    """
    Create a numpy array with the model data.
    """
    ref_genome = [i for i in list(aln_seqdict.keys()) if i.startswith(genome_id)]

    assert len(ref_genome) > 0, "{} must be in alignment file".format(genome_id)
    ref_genome = ref_genome[0]

    if genome_id not in set(aln_seqdict.keys()):
        warnings.warn(
            "{} is set as reference genome".format(ref_genome)
        )

    sequence = aln_seqdict[ref_genome].replace("-", "")

    # get features
    num_genes = get_gene_count(genbank)
    gc_content = GC(sequence)
    pos_change_var = get_position_variance(aln_seqdict)
    cai = calc_cai(sequence, genbank)

    # create and standardize model data array
    model_data = np.array([num_genes, gc_content, pos_change_var, cai])
    print(model_data)
    model_data = (model_data - model_params['scaler_mean']) / np.sqrt(model_params['scaler_var'])
    print(model_data)
    return model_data


def get_model_params(training_families=None):
    """
    Extract model parameters given families used for training
    """
    if training_families is None:
        model_params = MODEL_PARAMS['Coronaviridae']
    else:
        assert training_families in set(MODEL_PARAMS.keys()), \
            "{} is not a valid training_family variable".format(training_families)
        model_params = MODEL_PARAMS[training_families]

    return model_params


def calc_incubation(fasta_path, genbank_path, reference_genome, training_families='Coronaviridae'):
    """
    Applies the model to a user provided data
    """
    with open(fasta_path) as alignment_handle:
        aln_seqdict = proc_fasta(alignment_handle, alphabet="nucleotide")

    genbank = read_gb(genbank_path, reference_genome)
    model_params = get_model_params(training_families)

    model_data = create_model_data(
        reference_genome,
        aln_seqdict,
        genbank,
        model_params
    )

    ret_incubation = predict(model_data, model_params)
    return ret_incubation


# Main -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--alignment", type=str,
        help="Alignment of genomes in fasta format. Must contain an entry with the genome id.",
        required=True
    )

    parser.add_argument(
        "--genome-id", type=str, help="ID of the reference genome", required=True
    )

    parser.add_argument(
        "--genbank", type=str,
        help="Genbank file with the reference genome. Must contain an entry with the genome id.",
        required=True
    )

    parser.add_argument(
        "--training_families", type=str,
        help="Viral families used for model training. Options: "
             "A. 'Coronaviridae', "
             "B.'Coronaviridae_Paramyxoviridae', "
             "C. 'Coronaviridae_Pneumoviridae'. "
             "Default is 'Coronaviridae'",
        required=False
    )

    warnings.simplefilter('once', UserWarning)

    args = parser.parse_args()

    random.seed(3218)
    aln_seqdict = proc_fasta(open(args.alignment), alphabet="nucleotide")
    genbank = read_gb(args.genbank, args.genome_id)
    model_params = get_model_params(args.training_families)

    model_data = create_model_data(args.genome_id, aln_seqdict, genbank, model_params)

    print(args.genome_id, round(predict(model_data, model_params), 3))
