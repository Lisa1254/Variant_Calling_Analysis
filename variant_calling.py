#!/usr/env/bin/ python3

'''
Inputs expected are BAM file and variable_sites file.
If no index for BAM file exists in same folder as BAM file, this script will construct one
Outputs will return in working directory, which does not have to be the same as input file location
'''

import sys
import os
import argparse
import re
import pysam
import pandas as pd
import matplotlib.pyplot as plt

###
### Set up environment
###

argparser = argparse.ArgumentParser(
        description="Import and interpret bam file for loci of interest")
argparser.add_argument("-b", "--bamfile", action="store",
		help="Provide location of bam file for interpretation.")
argparser.add_argument("-l", "--locifile", action="store",
		help="Provide location of file describing loci of interest.")

args = argparser.parse_args(sys.argv[1:])

BAMFILE=args.bamfile
LOCIFILE=args.locifile

#Ensure that index is present
bam_filename = BAMFILE.split("/")[-1]
index_file = bam_filename + ".bai"

if not os.path.exists(index_file):
    pysam.index(BAMFILE)

# If desired, change name of output files
part_one_output = "reads_per_locus.txt"
part_two_output = "biallelic_locus_data.txt"

###
### Read in files
###

#Read file of variable sites
locifile_open = open(LOCIFILE)
locifile_lines = locifile_open.readlines()
locifile_open.close()

#Read bam file
samfile = pysam.AlignmentFile(BAMFILE, "rb")

###
### Code part one
###

#Function to parse CIGAR string for locus position in read
def cigar_read(start, loc, cigar):
    pos = 0
    dels = 0
    ins = 0
    cigar_exp = re.findall('(\d+)([IDSMN])', cigar)
    for tup in cigar_exp:
        if tup[1] in 'M':
            pos += int(tup[0])
            if (loc - start + 1) <= pos:
                locus_out = loc - start - dels + ins + 1
                return(locus_out)
        elif tup[1] in 'DN':
            dels += int(tup[0])
            pos += int(tup[0])
        elif tup[1] in 'IS':
            pos += int(tup[0])
            ins += int(tup[0])



#Use function to parse BAM file for required information and write to specified file
def reads_per_locus_write(locus, outfile):
    out_open = open(outfile, "a")
    locus_split = locus.split("\t")
    locus_chr = locus_split[0]
    locus_pos = int(locus_split[1])
    locus_ref = locus_split[2]
    for read in samfile.fetch(locus_chr, start=locus_pos-1, end=locus_pos):
        read_delim = str(read).split("\t")
        read_name = read_delim[0]
        read_len = len(read_delim[9])
        read_cigar_string = read_delim[5]
        start_pos = int(read_delim[3])
        if read_cigar_string == "None":
            locus_pos_in_read = locus_pos - start_pos + 1
        else:
            locus_pos_in_read = cigar_read(start_pos, locus_pos, read_cigar_string)
        read_base = read_delim[9][locus_pos_in_read -1]
        out_reads_per_locus_tuple = (locus_chr, str(locus_pos), locus_ref, read_name, str(read_len), read_cigar_string, str(locus_pos_in_read), read_base)
        out_reads_per_locus = "\t".join(out_reads_per_locus_tuple) + "\n"
        out_open.write(out_reads_per_locus)
    out_open.close()


#Initialize file with header
out_open = open(part_one_output, "w")
out_open.write("locus_chr\tlocus_pos\tlocus_ref\tread_name\tread_len\tread_cigar_string\tlocus_pos_in_read\tread_base\n")
out_open.close()

#Use function to append lines while iterating over all loci of interest
#Skip header line at locifile_lines[0]
for line in locifile_lines[1:]:
    test_locus = line.rstrip("\n")
    reads_per_locus_write(test_locus, part_one_output)




#At end of code section that utilizes the bam file, best practices from pysam
#recommends closing the AlignmentFile object
samfile.close()


###
### Code part two
###


#Import data from saved file
df_variants = pd.read_csv(part_one_output, sep="\t")

#Use function to parse dataframe of part 1 results for required information and write to file
#Inputs are locus of interest, dataframe of part 1 results grouped by locus and read base with counts, and outfile name
def biallelic_locus_write(locus, df_grouped, outfile):
    out_open = open(outfile, "a")
    locus_split = locus.split("\t")
    locus_pos = int(locus_split[1])
    sub = df_grouped.loc[df_grouped["locus_pos"] == locus_pos]
    if len(sub["read_base"]) == 2:
        var_count = sub["count"].sum()
        locus_chr = locus_split[0]
        ref_base = locus_split[2]
        base1 = sub.iloc[0,1]
        base1_ct = int(sub.iloc[0,2])
        base2 = sub.iloc[1,1]
        base2_ct = int(sub.iloc[1,2])
        if ref_base == base1:
            ref_freq = base1_ct / var_count
            alt_base = str(base2)
            alt_freq = base2_ct / var_count
        else:
            ref_freq = base2_ct / var_count
            alt_base = str(base1)
            alt_freq = base1_ct / var_count
        if (((ref_base in ("A", "G")) and (alt_base in ("A", "G"))) or ((ref_base in ("C", "T")) and (alt_base in ("C", "T")))):
            substitution_type = "transition"
        else:
            substitution_type = "transversion"
        output_biallelic_tuple = (locus_chr, str(locus_pos), ref_base, alt_base, str(ref_freq), str(alt_freq), substitution_type)
        output_biallelic = "\t".join(output_biallelic_tuple) + "\n"
        out_open.write(output_biallelic)
    out_open.close()



#Construct dataframe of part 1 results grouped by locus and read base with counts
df_gb_readbase = df_variants.groupby(["locus_pos", "read_base"], as_index=False)["read_base"].value_counts()

#Initialize outfile with header
out_open = open(part_two_output, "w")
out_open.write("locus_chr\tlocus_pos\tref_base\talt_base\tref_freq\talt_freq\tsubstitution_type\n")
out_open.close()

#Iterate over each locus
for line in locifile_lines[1:]:
    locus = line.rstrip("\n")
    biallelic_locus_write(locus, df_gb_readbase, part_two_output)


###
### Code part three
###

#Uses data from part one (imported part 2 as df_variants)
#Separate matching from non matching bases
df_match = df_variants.loc[df_variants["locus_ref"] == df_variants["read_base"]]
df_not_match = df_variants.loc[df_variants["locus_ref"] != df_variants["read_base"]]

#Extract data for histogram
match_pos = list(df_match.iloc[0:,6])
not_match_pos = list(df_not_match.iloc[0:,6])

#Save matching histogram
plt.hist(match_pos, bins=101, facecolor='blue', alpha=0.5, edgecolor='black', label="Match")
plt.ylabel('Counts')
plt.xlabel('Position in Read')
plt.legend()
plt.title("Read Base Matching Reference Base")

plt.savefig('position_base_matching_ref_maxbin.png')
plt.clf()

#Save not matching histogram
plt.hist(not_match_pos, bins=101, facecolor='red', alpha=0.5, edgecolor='black', label="Not Match")
plt.ylabel('Counts')
plt.xlabel('Position in Read')
plt.legend()
plt.title("Read Base Not Matching Reference Base")

plt.savefig('position_base_not_matching_ref_maxbin.png')
plt.clf()

#Save mixed histogram
plt.hist(match_pos, bins=101, facecolor='blue', alpha=0.5, edgecolor='black', label="Match")
plt.hist(not_match_pos, bins=101, facecolor='red', alpha=0.5, edgecolor='black', label="Not Match")
plt.ylabel('Counts')
plt.xlabel('Position in Read')
plt.legend()
plt.title("Comparison of Locus Position in Read")

plt.savefig('position_base_compare_maxbin.png')
plt.clf()

###
### Code part four
###

#This figure uses the data from part two
df_biallelic = pd.read_csv(part_two_output, sep="\t")

#Use pandas groupby to count substitution types
df_substitutions = df_biallelic.groupby(["substitution_type"], as_index=False)["substitution_type"].value_counts()

#Check which order groupby returns counts and save with "transitions" first
if df_substitutions.iloc[0,0] == "transition":
    counts = [int(df_substitutions.iloc[0,1])]
    counts.append(int(df_substitutions.iloc[1,1]))
else:
    counts = int(df_substitutions.iloc[1,1])
    counts.append(int(df_substitutions.iloc[0,1]))

#Define x axis labels
substitution = ['transition', 'transversion']

#Plot
plt.bar(substitution, counts, color ='maroon', width = 0.4)
plt.xlabel("Substitution Type")
plt.ylabel("Number of Base Substitutions")
plt.title("Substitution Type of Biallelic Reads")

#Save figure
plt.savefig('biallelic_substitution_types.png')
plt.clf()



