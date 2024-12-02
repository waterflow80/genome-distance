import Bio
from Bio import SeqIO
from Bio.UniProt.GOA import record_has

from config import GENOME_DATA_PATH, K_MERS_VAL
import os
from typing import Dict, List, Tuple
from collections import defaultdict
import hashlib
import itertools

def load_data()->Dict[str, str]:
    """
    Read the content of the fasta files and return a dict of sequence's ids mapped to the actual sequence
    :return:
    """
    genomes_dict = {}  # Map contains the sequence name of the file and the corresponding sequence
    os.chdir(GENOME_DATA_PATH)
    fasta_files = os.listdir()
    for fasta_file in fasta_files:
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        genomes_dict.update(record_dict)

    print(f"Returning fasta files {genomes_dict.keys()}")
    return genomes_dict

def apply_k_mers(genomes_dict:dict, k:int)->dict:
    """
    Split the genomes' sequences using the k-mers technique
    :param
    :param k:
    :return: {'seq_id1' : ["ATTCCGT", "ATTCAACT", ...],
              'seq_id2' : ["ATTCCGT", "ATTCAACT", ...],
            }
    """
    print(f"Applying k-mers..")
    genomes_k_mers = defaultdict(list)  # {seq_id: ["ATTGCC...TGCA", "ATTAAC...TGGC", ...]}
    for seq_id, sequence in genomes_dict.items():
        print(f"Calculating kmers for sequence {seq_id}")
        seqs = [] # the resulting sequences
        i = 0
        while (k + i) <= len(sequence):
            seqs.append(sequence.seq[i:k + i])
            i += 1
        genomes_k_mers[seq_id] = seqs

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

def calculate_jaccard_index(A: List, B: List) -> float:
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
    print(f"Calculating sketch of k-mers genomes..")
    for seq_id, sequences in k_mers_genomes.items():
        # step 1
        forward_hash = [calculate_hash(kmer) for kmer in sequences]
        reversed_compliment_hash = [calculate_hash(get_reverse_compliment(kmer)) for kmer in sequences]
        res_hash = [min(fwd_hash, cmp_rev_hash) for fwd_hash, cmp_rev_hash in zip(forward_hash, reversed_compliment_hash)]

        # step 2
        res_hash = sorted(res_hash)

        # step 3 update
        k_mers_genomes[seq_id] = res_hash[:1000]

    return k_mers_genomes

def get_all_sequences(k_mers_genomes):
    all_sequences = []  # List of tuples of sequences
    for genome, kmers in k_mers_genomes.items():
        for seq in kmers:
            all_sequences.append(seq)
    return all_sequences

if __name__ == '__main__':
    genomes = load_data()
    k_mers_genomes = apply_k_mers(genomes, k=K_MERS_VAL)
    k_mers_genomes_sketch = calculate_sketch(k_mers_genomes)

    #Calculating distance using the sketch
    all_sequences =list( k_mers_genomes_sketch.keys())
    for comb in itertools.combinations(all_sequences, 2):
        print(f"Distance between ", comb[0], " and ", comb[1], "is: ",
              calculate_jaccard_distance(k_mers_genomes_sketch.get(comb[0]), k_mers_genomes_sketch.get(comb[1])))

