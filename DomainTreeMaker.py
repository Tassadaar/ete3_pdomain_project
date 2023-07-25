import os
import sys
from Bio import SeqIO
from pyhmmer import easel, plan7, hmmer
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node
from Domain import Domain


def write_unaligned_sequences(sa):
    unaligned_seqs = sa.replace("-", "")

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


def generate_hits_file(sa):
    hits_file = "toy_topHits.txt"

    # unalign SA
    with open(sa) as sa_file:
        write_unaligned_sequences(str(sa_file.read()))

    alphabet = easel.Alphabet.amino()

    # loading the alignment
    with easel.SequenceFile("unaligned_toy.fasta", format="fasta", digital=True, alphabet=alphabet) as seq_file:
        sequences = seq_file.read_block()

    # load all the HMMs from an HMM file into a list
    with plan7.HMMFile("Pfam-A.hmm") as hmm_file:
        hmm = list(hmm_file)

    all_hits = list(hmmer.hmmscan(sequences, hmm, cpus=0))

    with open(hits_file, "wb") as output:
        for hits in all_hits:
            hits.write(output, "domains", True)

    # clean up temp files
    delete_file("unaligned_toy.fasta")

    return hits_file


def parse_hits_file(input_file):
    all_domains = {}

    with open(input_file, "r") as hits_file:

        current_parent = None
        for line in hits_file:

            if line[0] == "#":  # skipping headers
                continue

            words = line.split()
            domain = Domain(words[0], words[3], words[17], words[18])  # coordinates are assumed aligned
            if domain.parent != current_parent:
                all_domains[domain.parent] = [domain.get_motif_format()]
                current_parent = domain.parent
            else:
                all_domains[domain.parent].append(domain.get_motif_format())

    # delete_file("toy_topHits.txt")

    return all_domains


# for now, we are assuming the input file contains aligned amino acid sequences in fasta
def main():
    try:
        sa_file = "input_files/toy.fasta"
        tree_file = "input_files/toy.tree"

        hits_file = generate_hits_file(sa_file)

        domains = parse_hits_file(hits_file)

        tree = Tree(tree_file)
        leaves = tree.get_leaf_names()
        all_seqs = {}

        for seq_record in SeqIO.parse(sa_file, "fasta"):
            all_seqs[seq_record.id] = str(seq_record.seq)

        for leaf in leaves:
            seq_face = SeqMotifFace(all_seqs[leaf], seq_format="line", motifs=domains[leaf])
            (tree & leaf).add_face(seq_face, 0, "aligned")

        ts = TreeStyle()
        ts.show_scale = False
        tree.render("toy.png", h=50 * len(leaves), tree_style=ts)

    except FileNotFoundError:
        print(f"File not found!")


if __name__ == "__main__":
    main()
