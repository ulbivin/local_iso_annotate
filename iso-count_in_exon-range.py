import sys
import argparse
from collections import defaultdict

# Function to parse the annotated output file and sum occurrences within the exon range
def parse_and_sum_isoforms(file, exon_start, exon_end, strict):
    isoform_counts = defaultdict(int)
    
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            entry_name = parts[0]
            count = int(parts[1])
            exons = parts[2].split(', ')

            # Extract the range of exons specified
            range_exons = []
            my_start_ok = my_end_ok = False
            for exon in exons:
                if '-' in exon:
                    # Handle case where exon is in the form of 'start-end'
                    if my_start_ok and not my_end_ok:
                        range_exons.append(exon)
                else:
                    # Handle normal exon number
                    exon_num = int(exon.split('_')[0])
                    
                    if exon_start <= exon_num <= exon_end:
                        if strict:
                            range_exons.append(exon)
                        else:
                            range_exons.append(str(exon_num))
                    if exon_start == exon_num:
                        my_start_ok = True
                    if exon_end == exon_num:
                        my_end_ok = True

            # Only consider isoforms that have both exon_start and exon_end
            if my_start_ok and my_end_ok:
                range_exons_str = ', '.join(range_exons)
                isoform_counts[range_exons_str] += count


    return isoform_counts

# Main function to run the process
def main():
    parser = argparse.ArgumentParser(description='Sum isoform occurrences within a specified exon range.')
    parser.add_argument('-f', '--file', required=True, help='Annotated output file')
    parser.add_argument('-s', '--start', type=int, required=True, help='Start exon number')
    parser.add_argument('-e', '--end', type=int, required=True, help='End exon number')
    parser.add_argument('-S', '--strict', action='store_true', help='I set, full exon names will be used, if not, only exon numbers (default: False)')
    
    
    args = parser.parse_args()

    annotated_output_file = args.file
    exon_start = args.start
    exon_end = args.end
    optional_flag = args.strict

    isoform_counts = parse_and_sum_isoforms(annotated_output_file, exon_start, exon_end, optional_flag)

    for i, (exons, count) in enumerate(sorted(isoform_counts.items(), key=lambda item: -item[1]), 1):
        print(f"Local-Iso{i}\t{count}\t{exons}")

if __name__ == "__main__":
    main()
