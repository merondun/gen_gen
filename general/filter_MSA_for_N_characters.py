#!/usr/bin/env python3
#This script will take a multi sequence fasta and remove any sequences which contain more than -t "N" characters, e.g. script.py input.fa output.fa -t 0.2 removes sequences with >20% N's  
import argparse

def read_fasta(file):
    with open(file, 'r') as f:
        header, sequence = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, ''.join(sequence)
                header, sequence = line, []
            else:
                sequence.append(line)
        if header:
            yield header, ''.join(sequence)

def filter_sequences(input_file, output_file, threshold):
    with open(output_file, 'w') as output:
        for header, seq in read_fasta(input_file):
            len_seq = len(seq)
            n_count = seq.upper().count("N")
            n_percentage = n_count / len_seq
            if n_percentage <= threshold:
                output.write(header + "\n" + seq + "\n")

def main():
    parser = argparse.ArgumentParser(description="Filter sequences in a FASTA file based on the percentage of 'N' characters")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output FASTA file")
    parser.add_argument("-t", "--threshold", type=float, default=0.3, help="Threshold for 'N' character percentage (default: 0.3)")

    args = parser.parse_args()

    filter_sequences(args.input_file, args.output_file, args.threshold)

if __name__ == "__main__":
    main()


