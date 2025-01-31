# local_iso_annotate
************************************** UPDATE 20250131 **************************************

Updated with Readme and new exon info for BRCA1 annotation file.
Please note that iso-count_in_exon-range.py should no longer be used - use script by José Marcos instead.

Please also note updated ranges (systematic numbering):

| Mutation | Variant | Exons reange | Start Exon | End Exon |
|----------|---------|---------------|------------|----------|
| Mut1     | BRCA1 c.135-1G>T  | Exon 2 to Exon 6  | 2  | 6  |
| Mut2     | BRCA2 c.8632+1G>A | Exon 18 to Exon 21 | 18 | 21 |
| Mut3     | BRCA1 c.594-2A>C  | Exon 8 to Exon 12  | 7  | 13 |
| Mut4     | BRCA2 c.426-12_8delGTTTT | Exon 2 to Exon 8  | 2  | 8  |
| Mut5     | BRCA2 c.9501+3A>T | Exon 24 to Exon 27 | 22 | 26 |
| Mut6     | BRCA2 c.7988A>T   | Exon 16 to Exon 19 | 16 | 19 |
| Mut8     | BRCA1 c.671-2A>G  | Exon 8 to Exon 12  | 7  | 13 |
| Mut9     | BRCA1 c.5467+5G>C | Exon 22 to Exon 24 | 19 | 23 |


*************************************** End of update ***************************************


bed12_annotate.py takes in an exon annotation file and the [sample]_all_corrected.bed file from the FLAIR-output, and translate the bed-file entries to exon numbers if they match:
python bed12_annotate.py [exon_info.tsv] [input.bed12] > [annotated_output.txt] 
 
The script will allow the first and last exon entry to start and end anywhere within a given exon if it has the annotated donor or acceptor site.
 
comb_and_sort_annotation_counts.py combines the annotated output from the previous script with the [sample]_counts_matrix.tsv.counts.tsv file from the FLAIR-output:
python comb_and_sort_annotation_counts.py  [annotated_output.txt] [counts.tsv] > [annotated_counts.txt]
 
iso-count_in_exon-range.py is run on the annotated_counts.txt with a given range of exons and calculates the number of “local isoforms” within this:
python iso-count_in_exon-range.py -f [annotated_counts.txt] -s [start exon] -e [end exon]
 
The output is like this example looking at exon 4 to 10 in a BRCA2 sequencing:
python /mnt/DATA1/FLAIR_kConFab/iso-count_in_exon-range.py -f 61_071110075_BRCA2_ex1F-ex11R_seq040_counts_matrix.tsv.counts.tsv.anno_counts.tsv -s 4 -e 10
Local-Iso1      4380    4, 5, 6, 7, 8, 9, 10
Local-Iso2      913     4, 5, 8, 9, 10
Local-Iso3      503     4, 5, 6, 8, 9, 10
Local-Iso4      465     4, 5, 7, 8, 9, 10
Local-Iso5      356     4, 6, 7, 8, 9, 10
Local-Iso6      285     4, 8, 9, 10
Local-Iso7      278     4, 5, 6, 7, 32329443-32331030, 10
Local-Iso8      207     4, 7, 8, 9, 10
Local-Iso9      150     4, 5, 32329443-32331030, 10
 
The one/two digit numbers are identified exons, and the longer numbers denote genomic positions of observed exons not in the annotation file (chromosome number is omitted). 32329443-32331030 in the example corresponds to exon 8-start and exon 9-end meaning intron retention.
 
For BRCA1, I have added a number of alternative exons to the exon-info file e.g. the alternative 3nt shorter exon 13 is denoted 13_dE13p3 (d is for delta). More alternative exons can be added if desired – it is just important to have the number followed by underscore and no non-ASCII characters (like the delta sign) in the name.
 
iso-count_in_exon-range.py will treat the added alternative exons like 13_dE13p3 as exon 13, unless the strict [-S] parameter is set.
 
E.g. python /mnt/DATA1/FLAIR_kConFab/iso-count_in_exon-range.py -f 57_071110075_BRCA1_5UTR-3UTR_seq050_counts_matrix.tsv.counts.tsv.anno_counts.tsv -s 11 -e 16
Will print:
Local-Iso1      3739    11, 12, 13, 14, 15, 16
Local-Iso2      276     11, 12, 43079334-43079399, 13, 14, 15, 16
Local-Iso3      31      11, 13, 14, 15, 16
Local-Iso4      13      11, 12, 13, 15, 16
While python /mnt/DATA1/FLAIR_kConFab/iso-count_in_exon-range.py -f 57_071110075_BRCA1_5UTR-3UTR_seq050_counts_matrix.tsv.counts.tsv.anno_counts.tsv -s 11 -e 16 -S
The script differentiates between 13 and 13_Ex13p3:
Local-Iso1      2833    11, 12, 13, 14, 15, 16
Local-Iso2      906     11, 12, 13_Ex13p3, 14, 15, 16
Local-Iso3      203     11, 12, 43079334-43079399, 13, 14, 15, 16
Local-Iso4      73      11, 12, 43079334-43079399, 13_Ex13p3, 14, 15, 16
Local-Iso5      31      11, 13_Ex13p3, 14, 15, 16
Local-Iso6      13      11, 12, 13_Ex13p3, 15, 16

