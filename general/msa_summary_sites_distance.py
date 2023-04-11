import os
import sys
import glob
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator

def main(fasta_directory, ref_seq_id):
    fasta_files = glob.glob(os.path.join(fasta_directory, '*.fa'))

    if not fasta_files:
        print("No FASTA files found in the given directory.")
        return

    with open('output.tsv', 'w') as output_file:
        output_file.write('file_name\tsequence_id\tsequence_length\tnum_mutations\tnum_polymorphic\tnum_N_and_gaps\tbs_distance\n')

        for fasta_file in fasta_files:
            print(f"Processing {fasta_file}")
            ref_seq = None
            sequences = list(SeqIO.parse(fasta_file, 'fasta'))

            # Convert all sequences to uppercase
            for seq in sequences:
                seq.seq = seq.seq.upper()

            for seq in sequences:
                if seq.id == ref_seq_id:
                    ref_seq = seq
                    break

            if ref_seq is None:
                print(f"Reference sequence with ID {ref_seq_id} not found in {fasta_file}")
                continue

            polymorphic_sites = set()
            calculator = DistanceCalculator('blosum62')

            for seq in sequences:
                if seq.id != ref_seq_id:
                    for i, (a, b) in enumerate(zip(ref_seq, seq)):
                        if a != b and a not in {'N', '-'} and b not in {'N', '-'}:
                            polymorphic_sites.add(i)

            for seq in sequences:
                num_mutations = 0
                num_N_and_gaps = 0
                bs_distance = 0
                if seq.id != ref_seq_id:
                    num_mutations = sum(1 for a, b in zip(ref_seq, seq) if a != b and a not in {'N', '-'} and b not in {'N', '-'})
                    num_N_and_gaps = sum(1 for a, b in zip(ref_seq, seq) if (a == 'N' or a == '-') or (b == 'N' or b == '-'))
                    
                    # Create a MultipleSeqAlignment object for the pair of sequences and calculate the distance
                    pair_alignment = MultipleSeqAlignment([ref_seq, seq])
                    dm = calculator.get_distance(pair_alignment)
                    bs_distance = dm[ref_seq.id, seq.id]  

                output_file.write(f"{os.path.basename(fasta_file)}\t{seq.id}\t{len(seq)}\t{num_mutations}\t{len(polymorphic_sites)}\t{num_N_and_gaps}\t{bs_distance}\n")
                print(f"Writing mutation information for sequence {seq.id}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <fasta_directory> <reference_sequence_id>")
    else:
        main(sys.argv[1], sys.argv[2])

