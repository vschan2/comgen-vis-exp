import argparse
import csv
import os

from Bio import SeqIO


def process_genome(fasta_path, output_path):
    record = SeqIO.read(fasta_path, "fasta")

    filename = os.path.basename(fasta_path)

    # Create a CSV file
    csv_file = f"{output_path}/{filename}.csv"
    with open(csv_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        writer.writerow(["Position", "Base", "Count"])

        # Write the sequence data row by row
        for i, base in enumerate(record.seq, start=1):
            writer.writerow([i, base, 1])


def process_genome_2(fasta_path, output_path):
    record = SeqIO.read(fasta_path, "fasta")

    filename = os.path.basename(fasta_path)

    # test_read = [ "A", "T", "G", "C", "A", "A", "T", "C", "A", "A", "A" ]

    # Write the sequence data row by row
    index_dict = 0
    base_dict = {}
    for i, base in enumerate(record.seq, start=0):
    # for i, base in enumerate(test_read, start=0):
        if i == 0:
            prev_base = base
        
        if i == 0 or base != prev_base:
            base_dict[index_dict] = { "base": base, "start": i, "end": i+1, "count": 1 }
            index_dict += 1
        elif base == prev_base:
            base_dict[index_dict-1]["end"] += 1
            base_dict[index_dict-1]["count"] += 1
        
        prev_base = base
    
    # print(base_dict)

    # Create a CSV file
    csv_file = f"{output_path}/{filename}.csv"
    create_csv_file(csv_file, base_dict)


def process_genome_3(fasta_path, output_path):
    record = SeqIO.read(fasta_path, "fasta")

    filename = os.path.basename(fasta_path)

    # test_read = [ "A", "T", "G", "C", "A", "A", "T", "C", "A", "A", "A" ]

    # Write the sequence data row by row
    index_dict = 0
    base_dict = {}
    for i, base in enumerate(record.seq, start=0):
    # for i, base in enumerate(test_read, start=0):
        base_dict[i] = { "base": base, "start": i, "end": i+1, "count": 1 }
    
    # print(base_dict)

    # Create a CSV file
    csv_file = f"{output_path}/{filename}.csv"
    create_csv_file(csv_file, base_dict)


def create_csv_file(csv_file, base_dict):
    with open(csv_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        writer.writerow(["Base", "Start", "End", "Count"])

        for item_key, item_value in list(base_dict.items()):
            # print(item_value["base"])
            writer.writerow([item_value["base"], item_value["start"], item_value["end"], item_value["count"]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTA into CSV with each row represents the position and its base")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("--output", default=".", help="Output path of CSV file")
    
    args = parser.parse_args()
    
    process_genome_3(args.fasta_file, args.output)