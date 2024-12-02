import Bio
from Bio import SeqIO
from config import GENOME_DATA_PATH, K_MERS_VAL
import os
from typing import Dict, List
from collections import defaultdict
import hashlib
import itertools

def load_data()->Dict[str, SeqIO.FastaIO]:
    """
    Read the content of the fasta files and return a list of Bio.SeqIO.FastaIO.FastaIterator
    objects
    :return: dict of fasta_file_name -> SeqIO.FastaIO.FastaIterator
    """
    genomes_fasta_iterators_dict = {}  # Map contains the name of the file and the corresponding fastaIterator
    os.chdir(GENOME_DATA_PATH)
    fasta_files = os.listdir()
    for fasta_file in fasta_files:
        fasta_iter = SeqIO.parse(fasta_file, "fasta")
        genomes_fasta_iterators_dict[fasta_file] = fasta_iter

    print(f"Returning fasta files {genomes_fasta_iterators_dict.keys()}")
    return genomes_fasta_iterators_dict


def get_full_sequence_from_fasta_iterator(fasta_iterator:SeqIO.FastaIO.FastaIterator)->str:
    """
    Concatenate all the sequences in the given fasta file and return one sequence
    :param fasta_iterator:
    :return:
    """
    print(f"Extracting full sequence from {fasta_iterator} ..")
    full_sequence = ""
    for record in list(fasta_iterator):
        full_sequence += record.seq
    return full_sequence


def apply_k_mers(genomes_fasta_iterators:Dict[str, SeqIO.FastaIO], k:int):
    """
    Split the genomes' sequences using the k-mers technique
    :param genomes_fasta_iterators:
    :param k:
    :return:
    """
    print(f"Applying k-mers..")
    genomes_k_mers = defaultdict(list)  # {fasta_filename: ["ATTGCC...TGCA", "ATTGCC...TGCA", ...]}
    for fasta_file, fasta_iterator in genomes_fasta_iterators.items():
        full_sequence = get_full_sequence_from_fasta_iterator(fasta_iterator)
        i = 0
        while (k + i) <= len(full_sequence):
            genomes_k_mers[fasta_file].append(full_sequence[i:k + i])
            i += 1
    return genomes_k_mers

def _get_intersection(A, B)->set:
    """
    Return a set of intersection of A and B
    :param A:
    :param B:
    :return:
    """
    return set(A).intersection(set(B))

def _get_union(A, B):
    """
    Calculate the union of A and B
    :param A:
    :param B:
    :return:
    """
    return set.union(set(A), set(B))

def calculate_jaccard_index(A: List[Bio.Seq.Seq], B: List[Bio.Seq.Seq]) -> float:
    """
    Calculate the Jaccard distance between A and B
    :param A:
    :param B:
    :return:
    """
    return len(_get_intersection(A, B)) / len(_get_union(A, B))

def calculate_jaccard_distance(A, B):
    return 1 - calculate_jaccard_index(A, B)

def calculate_hash(seq, seed=0):
    return hashlib.md5(str(seq).encode()).hexdigest()

def get_reverse_compliment(seq):
    """
    Return the reverse compliment of the given sequence
    :param seq:
    :return:
    """
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}  # Some of the sequences contains 'N' (unknown/uncertain)
    reverse_complement = "".join(complements[base] for base in reversed(str(seq)))
    return reverse_complement

def calculate_sketch(k_mers_genomes:dict):
    """
    Calculate the sketch of the given k-mers genomes. The sketch is obtained following these steps:
    1. Converting all 14-mers from the genome into integers using a hash function. For each individual k-mer, a
    hash is generated for its forward and reverse complement, and the lowest value out of these two is chosen as
    the canonical k-mer.
    2. Sorting these hashes.
    3. Taking the smallest S hashes, where s=1000 is the ‘sketch size’. This can be made larger or smaller,
    but for this purpose we’ll use a fixed value of 1000.
    :param k_mers_genomes:
    :return:
    """
    for genome, kmers in k_mers_genomes.items():
        # step 1
        forward_hash = [calculate_hash(kmer) for kmer in kmers]
        reversed_compliment_hash = [calculate_hash(get_reverse_compliment(kmer)) for kmer in kmers]
        res_hash = [min(fwd_hash, cmp_rev_hash) for fwd_hash, cmp_rev_hash in zip(forward_hash, reversed_compliment_hash)]

        # step 2
        res_hash = sorted(res_hash)

        # step 3
        k_mers_genomes[genome] = res_hash[:1000]

    return k_mers_genomes

if __name__ == '__main__':
    genomes = load_data()
    k_mers_genomes = apply_k_mers(genomes, k=K_MERS_VAL)

    #Calculating distance using the sketch
    fasta_files = list(k_mers_genomes.keys())
    for comb in itertools.combinations(fasta_files, 2):
        print(f"Distance between ", comb[0], " and ", comb[1], "is: ",
              calculate_jaccard_distance(k_mers_genomes.get(comb[0]), k_mers_genomes.get(comb[1])))

