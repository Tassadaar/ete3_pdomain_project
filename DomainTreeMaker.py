import sys
import pyhmmer
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node


# for now, we are assuming the input file contains amino acid sequences in fasta
def main():
    alphabet = pyhmmer.easel.Alphabet.amino()

    # loading the alignment
    with pyhmmer.easel.MSAFile("input_files/nematode.fasta", format="afa", digital=True, alphabet=alphabet) as msa_file:
        msa = list(msa_file.read().sequences)

    # load all the HMMs from an HMM file into a list
    with pyhmmer.plan7.HMMFile("Pfam-A.hmm") as hmm_file:
        hmm = list(hmm_file)

    all_hits = list(pyhmmer.hmmer.hmmscan(msa, hmm, cpus=0))

    with open("toy_topHits.txt", "wb") as output:
        for hits in all_hits:
            hits.write(output, "domains", True)


if __name__ == "__main__":
    main()
