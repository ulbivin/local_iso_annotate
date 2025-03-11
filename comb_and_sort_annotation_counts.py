import pandas as pd
import sys

# Function to parse the annotated output file
def parse_annotated_output(file):
    annotations = {}
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            entry_name = parts[0].rsplit('_' ,1)[0]
            exons = parts[1].strip().split(', ')
            annotations[entry_name] = exons
    return annotations

# Function to parse the counts TSV file
def parse_counts_tsv(file):
    counts = pd.read_csv(file, sep='\t', skiprows=1, header=None, names=["entry", "count"])
    counts["entry"] = counts["entry"].apply(lambda x: x.rsplit('_', 1)[0])
    counts_dict = dict(zip(counts["entry"], counts["count"].astype(int)))
    return counts_dict

# Function to combine and sort annotations with counts
def combine_and_sort_annotations(annotations, counts):
    combined_list = []
    for entry, exons in annotations.items():
        if entry in counts:
            count = counts[entry]
            exons_str = ', '.join(exons)
            combined_list.append((entry, count, exons_str))
    
    # Sort by count in descending order
    combined_list.sort(key=lambda x: x[1], reverse=True)
    
    return combined_list

# Main function to run the process
def main():
    if len(sys.argv) != 3:
        print("Usage: python combine_and_sort_annotations.py <annotated_output.txt> <counts.tsv>")
        sys.exit(1)
    
    annotated_output_file = sys.argv[1]
    counts_file = sys.argv[2]
    
    annotations = parse_annotated_output(annotated_output_file)
    counts = parse_counts_tsv(counts_file)
    combined_list = combine_and_sort_annotations(annotations, counts)
    
    for entry, count, exons_str in combined_list:
        print(f"{entry}\t{count}\t{exons_str}")

if __name__ == "__main__":
    main()
