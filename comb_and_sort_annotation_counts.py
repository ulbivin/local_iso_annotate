import pandas as pd
import sys

# Function to parse the counts file
def parse_counts(counts_file):
    """
    Reads the counts file into a pandas DataFrame.

    Parameters:
        counts_file (str): Path to the counts file.

    Returns:
        DataFrame: Parsed counts data.
    """
    counts = pd.read_csv(counts_file, sep="\t")
    return counts

# Function to parse the GTF file and extract exon information
def parse_gtf(gtf_file):
    """
    Reads a GTF file, extracts exon features, and removes 'chr' prefix from chromosome numbers.

    Parameters:
        gtf_file (str): Path to the GTF file.

    Returns:
        DataFrame: Filtered GTF data containing only exons.
    """
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None)
    gtf.columns = ["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    
    # Filter only exon entries
    gtf = gtf[gtf["feature"] == "exon"]
    
    # Remove 'chr' prefix from chromosome names
    gtf["chromosome"] = gtf["chromosome"].astype(str).str.replace("chr", "")
    
    return gtf

# Function to parse the exon information file
def parse_exon_info(exon_info_file):
    """
    Reads the exon information file and removes 'chr' from chromosome names.

    Parameters:
        exon_info_file (str): Path to the exon information file.

    Returns:
        DataFrame: Parsed exon information.
    """
    exon_info = pd.read_csv(exon_info_file, sep="\t")
    
    # Remove 'chr' prefix from chromosome names
    exon_info["chromosome"] = exon_info["chromosome"].astype(str).str.replace("chr", "")
    
    return exon_info

# Function to parse the BED12 annotation file
def parse_bed12_annotation(bed12_annotation_file):
    """
    Reads the BED12 annotation file into a pandas DataFrame.

    Parameters:
        bed12_annotation_file (str): Path to the BED12 annotation file.

    Returns:
        DataFrame: Parsed BED12 annotation data.
    """
    bed12_annot = pd.read_csv(bed12_annotation_file, sep="\t", header=None, names=["id", "exons"])
    return bed12_annot

# Function to map genomic coordinates to exon numbers
def map_exon_numbers(exon_coordinates, exon_info, chromosome):
    """
    Replaces genomic coordinates with exon numbers if an exact match is found.
    If no match is found, retains 'chromosome:start-end' format.

    Parameters:
        exon_coordinates (list of tuples): List of (start, end) coordinates.
        exon_info (DataFrame): DataFrame containing exon information.
        chromosome (str): Chromosome number.

    Returns:
        str: Comma-separated list of exon numbers or genomic coordinates.
    """
    mapped_exons = []
    
    for start, end in exon_coordinates:
        # Find an exact match in exon_info
        match = exon_info[(exon_info["start"] == start) & (exon_info["end"] == end)]
        
        if not match.empty:
            exon_number = match["exon#"].values[0]
            mapped_exons.append(str(exon_number))  # Append the exon number if a match is found
        else:
            mapped_exons.append(f"{chromosome}:{start}-{end}")  # Keep genomic coordinates if no match
    
    return ", ".join(mapped_exons)

# Function to annotate the counts file
def annotate_counts(counts, bed12_annot, exon_info, gtf_info):
    """
    Annotates the counts file with exon or transcript information.

    Parameters:
        counts (DataFrame): Parsed counts data.
        bed12_annot (DataFrame): Parsed BED12 annotation data.
        exon_info (DataFrame): Parsed exon information.
        gtf_info (DataFrame): Parsed GTF exon data.

    Returns:
        List[str]: Sorted annotated entries in formatted string representation.
    """
    annotated_entries = []
    
    for _, row in counts.iterrows():
        entry_id = row.iloc[0].strip()  # Clean up entry ID
        count = int(row.iloc[1])  # Extract counts
        counts_entry = entry_id.rsplit("_", 1)[0].strip().replace('_', '-')  # Format entry name
        
        # If the ID contains "ENST", it's a transcript
        if "ENST" in entry_id:
            transcript_id = entry_id.split("_")[0]
            if "-" in transcript_id:
                transcript_id = transcript_id.split("-")[0]
            # Find matching exons in GTF file
            matching_exons = gtf_info[gtf_info["attributes"].str.contains(f'transcript_id "{transcript_id}"', regex=False)]
            
            if not matching_exons.empty:
                # Extract exon start and end coordinates
                exon_coordinates = [(row["start"], row["end"]) for _, row in matching_exons.iterrows()]
                chromosome = matching_exons["chromosome"].values[0]  # Extract chromosome number
                
                # Map genomic coordinates to exon numbers
                exon_list = map_exon_numbers(exon_coordinates, exon_info, chromosome)
            else:
                exon_list = "N/A"  # If no match is found, mark as "N/A"
                sys.stderr.write(transcript_id+" not found in gtf file\n")
        
        else:
            # Find a match in BED12 annotation
            if len(counts_entry.rsplit("-", 1)[-1])==1:
                counts_entry=counts_entry.rsplit("-", 1)[0]
            match = bed12_annot[bed12_annot["id"].str.strip() == counts_entry]
            if not match.empty:
                exon_list = match["exons"].values[0]  # Get exon list from BED12 annotation
            else:
                exon_list = "N/A"  # If no match is found, mark as "N/A"
                sys.stderr.write(counts_entry+" not found in annotated bed12 tsv file\n")
        
        # Store annotated entry
        annotated_entries.append((counts_entry, count, exon_list))
    
    # Sort entries by counts (descending order)
    annotated_entries.sort(key=lambda x: x[1], reverse=True)

    # Format each entry as a tab-separated string
    return [f"{entry[0]}\t{entry[1]}\t{entry[2]}" for entry in annotated_entries]

# Main function to execute the script
def main():
    """
    Main function to parse input files, annotate counts, and print the results.
    """
    # Ensure the correct number of arguments are provided
    if len(sys.argv) != 5:
        print("Usage: python comb_and_sort_annotation_counts.py <counts.tsv> <bed12_annotated.tsv> <exon_info.tsv> <input.gtf>")
        sys.exit(1)
    
    counts_file = sys.argv[1]  # Path to counts file
    bed12_annotation_file = sys.argv[2]  # Path to BED12 annotation file
    exon_info_file = sys.argv[3]  # Path to exon information file
    gtf_file = sys.argv[4]  # Path to GTF file
    
    # Parse input files
    counts = parse_counts(counts_file)
    bed12_annot = parse_bed12_annotation(bed12_annotation_file)    
    exon_info = parse_exon_info(exon_info_file)
    gtf_info = parse_gtf(gtf_file)
    
    # Annotate counts
    annotated_data = annotate_counts(counts, bed12_annot, exon_info, gtf_info)
    
    with open(sys.argv[1]+'.anno_counts.tsv', 'w') as fout:
    # Print the annotated data
        for annotation in annotated_data:
            print(annotation, file=fout)

# Run main function if script is executed directly
if __name__ == "__main__":
    main()
