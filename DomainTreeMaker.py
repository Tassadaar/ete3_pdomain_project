import os
import sys
from pyhmmer import easel, plan7, hmmer
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node
from Domain import Domain


def write_unaligned_sequences(msa):
    unaligned_seqs = msa.replace("-", "")

    with open("unaligned_toy.fasta", "w") as out_file:
        out_file.write(unaligned_seqs)


def delete_file(filename):
    try:
        os.remove(filename)
        print(f"File '{filename}' has been deleted as part of clean-up.")
    except FileNotFoundError:
        print(f"File '{filename}' not found, so it could not be deleted.")
    except Exception as e:
        print(f"An error occurred while trying to delete the file '{filename}': {e}")


def parse_hits_file(input_file):
    domains = []

    with open(input_file, "r") as hits_file:
        try:

            for line in hits_file:

                if line[0] == "#":  # skipping headers
                    continue

                words = line.split()
                domains.append(Domain(words[0], words[3], words[17], words[18]))  # coordinates are assumed aligned

        except FileNotFoundError:
            print(f"File {input_file} not found!")

    return domains


# for now, we are assuming the input file contains aligned amino acid sequences in fasta
def main():
    # unalign MSA
    with open("input_files/nematode.fasta") as msa_file:
        write_unaligned_sequences(str(msa_file.read()))

    alphabet = easel.Alphabet.amino()

    # loading the alignment
    with easel.SequenceFile("unaligned_toy.fasta", format="fasta", digital=True, alphabet=alphabet) as seq_file:
        sequences = seq_file.read_block()

    # load all the HMMs from an HMM file into a list
    with plan7.HMMFile("Pfam-A.hmm") as hmm_file:
        hmm = list(hmm_file)

    all_hits = list(hmmer.hmmscan(sequences, hmm, cpus=0))

    with open("toy_topHits.txt", "wb") as output:
        for hits in all_hits:
            hits.write(output, "domains", True)

    # clean up temp files
    delete_file("unaligned_toy.fasta")


if __name__ == "__main__":
    main()
