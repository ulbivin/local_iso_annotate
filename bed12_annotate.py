import pandas as pd
import sys

# Function to parse the exon information file
def parse_exon_info(exon_info_file):
    """
    Reads the exon_info file into a pandas DataFrame.

    Parameters:
        exon_info_file (str): Path to the exon information file.

    Returns:
        DataFrame: Parsed exon information containing chromosome, start, end, exon name, exon number, and strand.
    """
    exon_info = pd.read_csv(exon_info_file, sep="\t")
    return exon_info

# Function to parse the BED12 file
def parse_bed12(bed12_file):
    """
    Reads the BED12 file into a pandas DataFrame.

    Parameters:
        bed12_file (str): Path to the BED12 file.

    Returns:
        DataFrame: Parsed BED12 data with standard BED12 format columns.
    """
    bed12 = pd.read_csv(bed12_file, sep="\t", header=None, names=[
        "chrom", "chromStart", "chromEnd", "name", "score", "strand", 
        "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
    ])
    return bed12

# Function to annotate BED12 entries with exon information
def annotate_bed12(exon_info, bed12):
    """
    Annotates BED12 entries by mapping block coordinates to known exons.

    Parameters:
        exon_info (DataFrame): DataFrame containing exon information.
        bed12 (DataFrame): DataFrame containing parsed BED12 data.

    Returns:
        List[str]: Annotated BED12 entries in a formatted string representation.
    """
    annotations = []
    
    # Iterate over each row in the BED12 file
    for index, row in bed12.iterrows():
        entry_name = row['name'].replace("_", "-")  # Replace underscores with hyphens in entry names
        chrom = row['chrom'].replace("chr", "")  # Remove "chr" prefix from chromosome
        strand = row['strand']
        chrom_start = row['chromStart']
        chrom_end = row['chromEnd']
        
        # Parse block starts and sizes (representing exon coordinates within the BED entry)
        block_starts = list(map(int, row['blockStarts'].strip(',').split(',')))
        block_sizes = list(map(int, row['blockSizes'].strip(',').split(',')))
        
        exon_list = []
        
        # Iterate over each block (exon) in the BED12 entry
        for i, (block_start, block_size) in enumerate(zip(block_starts, block_sizes)):
            block_start_pos = chrom_start + block_start + 1  # Convert 0-based BED start to 1-based
            block_end_pos = block_start_pos + block_size - 1  # Calculate block end position
            
            block_annotated = my_first = my_last = False  # Flags for annotation tracking
            
            # Iterate over the exon information to find a match
            for _, exon in exon_info.iterrows():
                exon_chrom = str(exon['chromosome'])  # Get chromosome from exon info
                exon_start = exon['start']  # Get exon start position
                exon_end = exon['end']  # Get exon end position
                exon_number = exon['exon#']  # Get exon number
                exon_name = exon['exon_name']  # Get exon name
                
                # Check if exon belongs to the same chromosome
                if exon_chrom == chrom:
                    # Check for first block match
                    if i == 0 and not my_first:
                        if exon_start <= block_start_pos and exon_end == block_end_pos:
                            exon_list.append(str(exon_number))  # Append exon number
                            block_annotated = my_first = True  # Mark first exon as annotated
                    # Check for last block match
                    elif i == len(block_starts) - 1 and not my_last:
                        if exon_start == block_start_pos and exon_end <= block_end_pos:
                            exon_list.append(str(exon_number))
                            block_annotated = my_last = True  # Mark last exon as annotated
                    # Check for middle block match
                    else:
                        if exon_start == block_start_pos and exon_end == block_end_pos:
                            exon_list.append(str(exon_number))
                            block_annotated = True
            
            # If no matching exon is found, store genomic coordinates instead
            if not block_annotated:
                exon_list.append(f"{chrom}:{block_start_pos}-{block_end_pos}")

        # Reverse the exon list if the strand is negative
        if strand == '-':
            exon_list.reverse()
        
        # Format annotation as "entry_name    exon1, exon2, exon3..."
        annotation = f"{entry_name}\t" + ", ".join(exon_list)
        annotations.append(annotation)
    
    return annotations

# Main function to execute the script
def main():
    """
    Main script execution function.

    Reads command-line arguments, parses the input files, annotates the BED12 entries,
    and prints the results.
    """
    # Ensure the correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python annotate_bed12.py <exon_info.tsv> <input.bed12>")
        sys.exit(1)
    
    exon_info_file = sys.argv[1]  # Path to exon info file
    bed12_file = sys.argv[2]  # Path to BED12 file
    
    # Parse input files
    exon_info = parse_exon_info(exon_info_file)
    bed12 = parse_bed12(bed12_file)
    
    # Annotate BED12 entries
    annotations = annotate_bed12(exon_info, bed12)
    
    # Print annotated entries
    for annotation in annotations:
        print(annotation)

# Run main function if script is executed directly
if __name__ == "__main__":
    main()
