import pandas as pd
import sys

# Function to parse the exon information file
def parse_exon_info(exon_info_file):
    exon_info = pd.read_csv(exon_info_file, sep="\t")
    return exon_info

# Function to parse the BED12 file
def parse_bed12(bed12_file):
    bed12 = pd.read_csv(bed12_file, sep="\t", header=None, names=[
        "chrom", "chromStart", "chromEnd", "name", "score", "strand", 
        "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
    ])
    return bed12

# Function to annotate the BED12 file with exon information
def annotate_bed12(exon_info, bed12):
    annotations = []
    
    for index, row in bed12.iterrows():
        entry_name = row['name']
        chrom = row['chrom'][3:]
        strand = row['strand']
        chrom_start = row['chromStart']
        chrom_end = row['chromEnd']
        block_starts = list(map(int, row['blockStarts'].strip(',').split(',')))
        block_sizes = list(map(int, row['blockSizes'].strip(',').split(',')))
        
        exon_list = []
        for i, (block_start, block_size) in enumerate(zip(block_starts, block_sizes)):
            block_start_pos = chrom_start + block_start + 1
            block_end_pos = block_start_pos + block_size - 1
            
            block_annotated = my_first = my_last = False
            
            for _, exon in exon_info.iterrows():
                exon_chrom = str(exon['chromosome'])
                exon_start = exon['start']
                exon_end = exon['end']
                exon_number = exon['exon#']
                exon_name = exon['exon_name']
                
                if exon_chrom == chrom:
                    if i == 0 and my_first == False:
                        # First: include overlapping exons
                        if exon_start <= block_start_pos and exon_end == block_end_pos:
                            exon_list.append(str(exon_number))
                            block_annotated = my_first = True
                    elif i == len(block_starts) - 1 and my_last == False:
                        # Last block: include overlapping exons
                        if exon_start == block_start_pos and exon_end <= block_end_pos:
                            #print (exon_start, block_start_pos, exon_end, block_end_pos)
                            exon_list.append(str(exon_number))
                            block_annotated = my_last = True
                    else:
                        # Intermediate blocks: only include exact exon boundaries
                        if exon_start == block_start_pos and exon_end == block_end_pos:
                            exon_list.append(str(exon_number))
                            block_annotated = True
            
            if not block_annotated:
                exon_list.append(str(block_start_pos)+'-'+str(block_end_pos))
        
        if strand == '-':
            exon_list.reverse()
        
        annotation = f"{entry_name}\t" + ", ".join(exon_list)
        annotations.append(annotation)
    
    return annotations

# Main function to run the annotation process
def main():
    if len(sys.argv) != 3:
        print("Usage: python annotate_bed12.py <exon_info.tsv> <input.bed12>")
        sys.exit(1)
    
    exon_info_file = sys.argv[1]
    bed12_file = sys.argv[2]
    
    exon_info = parse_exon_info(exon_info_file)
    bed12 = parse_bed12(bed12_file)
    annotations = annotate_bed12(exon_info, bed12)
    
    for annotation in annotations:
        print(annotation)

if __name__ == "__main__":
    main()
