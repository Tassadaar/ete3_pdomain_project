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


def generate_hits(sa):
    # unalign SA
    with open(sa) as sa_file:
        write_unaligned_sequences(str(sa_file.read()))

    alphabet = easel.Alphabet.amino()

    # loading the alignment
    with easel.SequenceFile("unaligned_toy.fasta", format="fasta", digital=True, alphabet=alphabet) as seq_file:
        sequences = seq_file.read_block()

    # load all the HMMs from an HMM file into a list
    with plan7.HMMFile("Pfam-A.h3m") as hmm_file:
        hmm = list(hmm_file)

    return hmmer.hmmscan(sequences, hmm, cpus=0)


def parse_hits_file(input_file):
    all_domains = {}

    with open(input_file, "r") as hits_file:

        current_parent = None
        for line in hits_file:

            if line[0] == "#":  # skipping headers
                continue

            words = line.split()
            domain = Domain(words[0], words[3], words[13], words[17], words[18])  # coordinates are assumed aligned

            if float(domain.score) < 50:  # skipping poor hits
                continue

            if domain.parent != current_parent:
                all_domains[domain.parent] = [domain.get_motif_format()]
                current_parent = domain.parent
            else:
                all_domains[domain.parent].append(domain.get_motif_format())

    # delete_file("toy_topHits.txt")

    return all_domains


def generate_domain_layout(start, end, name):

    if name not in static_colour_dict.keys():
        static_colour_dict[name] = colour_schemes.pop()

    return [start, end, "()", 100, 10, "black", static_colour_dict[name], f"arial|7|white|{name}"]


# for now, we are assuming the input file contains aligned amino acid sequences in fasta
def main():
    try:
        sa_file = "input_files/toy.fasta"
        tree_file = "input_files/toy.tree"

        all_hits = generate_hits(sa_file)

        domains = {}

        for hits in all_hits:
            for hit in hits:
                for doamin in hit.domains:
                    if doamin.score < 50:
                        continue

                    if hits.query_name.decode() not in domains.keys():
                        domains[hits.query_name.decode()] = [generate_domain_layout(doamin.env_from,
                                                                                    doamin.env_to,
                                                                                    hit.name.decode())]
                    else:
                        domains[hits.query_name.decode()].append(generate_domain_layout(doamin.env_from,
                                                                                        doamin.env_to,
                                                                                        hit.name.decode()))

        tree = Tree(tree_file)
        leaves = tree.get_leaf_names()
        all_seqs = {}

        for seq_record in SeqIO.parse(sa_file, "fasta"):
            all_seqs[seq_record.id] = str(seq_record.seq)

        for leaf in leaves:

            if leaf not in domains.keys():  # skipping sequences with no hits
                continue

            seq_face = SeqMotifFace(all_seqs[leaf], seq_format="line", motifs=domains[leaf])
            (tree & leaf).add_face(seq_face, 0, "aligned")

        tree_style = TreeStyle()
        tree_style.show_scale = False

        delete_file("unaligned_toy.fasta")

        tree.render("toy.png", h=50 * len(leaves), tree_style=tree_style)

    except FileNotFoundError:
        print(f"File not found!")


if __name__ == "__main__":
    colour_schemes = ["rgradient:blue", "rgradient:red", "rgradient:green", "rgradient:orange", "rgradient:purple"]
    static_colour_dict = {}
    main()
