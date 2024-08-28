#!/usr/bin/env python
# coding: utf-8

import os
#form data types in file. SeqIO handles sequence input/output
from Bio import SeqIO
from Bio.Seq import Seq
#import all sequence records. each gene file is read as one record (Locus to end of Origin)
from Bio.SeqRecord import SeqRecord
#feature location = datatype. categorize and allow parsing through features section
from Bio.SeqFeature import SeqFeature, FeatureLocation
#regular expressions
import re
#csv
import csv

#read record by record (one gene file at a time)
recs = [rec for rec in SeqIO.parse("genbankfiles_justkofirst.gb", "genbank")]


#number of records
print("Number of records is: ", len(recs))

#list for records that have target exons
result = []
#list for
exon_locations = []
#empty list for potential target exons that don't fit within critical region of particular gene
bad_exon = []
count = 0
#list for gene names that don't have 'target exon' in their record (untargeted genes)
untargeted_genes = []

#for loop going through each record
for rec in recs:
    #separate misc features into datatype. allow parsing of only misc features for start:end values
    feats = [feat for feat in rec.features if feat.type == "misc_feature"]
    #separate exon into datatype. allow parsing of only data under term 'exon'
    exons = [feat for feat in rec.features if feat.type == "exon"]
    # formatting feature key value
    for feat in feats:
      #look for mention of 'critical region' in misc_features datatype. always in 'note' section
      if(feat.qualifiers['note'][0].lower() == "critical region"):
        #print("Start :", feat.location.start)
        crit_start = feat.location.start
        #print("End :", feat.location.end)
        crit_end = feat.location.end

    crit_exon = []
    for exon in exons:

        #test for ccdc50
#        if (rec.description.split(":")[0].split(" ")[-1] == "Tm4sf4"):
#            print(exon.qualifiers)

        if("target exon" in exon.qualifiers['note'][0].lower()):
        #------------------------------------------
            crit_exon.append(exon.qualifiers['note'][0].split(' ')[3])
            if not (exon.location.start > crit_start and exon.location.end < crit_end):
              #print("Exon :", exon)
              bad_exon.append(str(exon.qualifiers['note'][0].split(' ')[3]))
              #target_exon.append(exon.qualifiers['note'][0].split(' ')[3])
              #result.append(exon.qualifiers['note'][0].split(' ')[3])

    #logging records that don't have target exon.
    if (len(crit_exon) == 0):
        #search through 'description' in the file and split by commas. append gene name (4th split) into new file
        untargeted_genes.append(rec.description.split(",")[2].split(" ")[4])
        continue

    #checking for multiple target exons within 1 gene. taking only first and last exon ID
    elif (len(crit_exon) > 1):
      result.append(crit_exon[0])
      exon_locations.append([crit_exon[0],crit_exon[-1]])
      count += 2
    #if gene has only 1 'target exon', append the ensembl ID
    else:
      result.append(crit_exon[0])
      exon_locations.append(crit_exon[0])
      count += 1
#count is how many mentions of 'target exon' there are. takes into consideration genes that have multiple target exons

#print("total mention of 'target exons':", count)
print("targeted genes: ", len(exon_locations))
print("untargeted genes: ", len(untargeted_genes))

#check
if 'ENSMUSE00000155698' in result:
    print('1')

#Writing all critical exon IDs into a file
file = open("crit_exon.txt","w")
for i in range(len(exon_locations)):
  file.write("{} = {}\n".format(result[i],exon_locations[i]))
file.close()

#Writing critical exon IDs that don't fit in critical range into a file
file = open("crit_exon_range.txt","w")
for i in bad_exon:
    file.write(i+"\n")
file.close()

#Writing untargeted gene names into a file
file = open("untargeted.txt","w")
for i in untargeted_genes:
    file.write(i+"\n")
file.close()


#read in BioMart file containing canonical data on all protein coding genes.
#format:Gene stable ID,Transcript stable ID,Exon stable ID,Exon region start (bp),Exon region end (bp),Start phase,End phase,Exon rank in transcript,Gene name,Chromosome/scaffold name
exons = open('pcg_canonical.txt','r').readlines()
#dictionary for gene IDs, gene names, ranks, trans IDs
gene_ids = {}
gene_names = {}
rank = {}
trans_ids_exon_first = {}
from collections import defaultdict

exon_dict = {}
for i in range(1,len(exons)):
  temp = exons[i].split(',')
  exon_dict.setdefault(temp[2], [])
  exon_dict[temp[2]] = (temp[1],int(temp[5]),int(temp[6]))
  gene_ids[temp[2]] = temp[0]
  trans_ids_exon_first[temp[2]] = temp[1]
  gene_names[temp[2]] = temp[8]
  rank[temp[2]] = int(temp[7])


#list for critical exon IDs from IMPC file that are not present in BioMart pcg_canonical file
non_canonical_exons = []

#fixing ordering of ids. using pcg_canonical file to order the genes with multiple critical exons in correct order (based on rank)
#need correct order to get right start and end phases
for i,exonid in enumerate(result):
    try:
        if (type(exon_locations[i]) == list):
            if (rank[exon_locations[i][0]] > rank[exon_locations[i][1]]):
                exon_locations[i][0], exon_locations[i][1] = exon_locations[i][1], exon_locations[i][0]
                result[i] = exon_locations[i][0]

    #append exon IDs that are not present in pcg_canonical into a new file
    except KeyError:
        non_canonical_exons.append(exonid)

print("Number of targeted exons not present in Ensembl canonical file: ", len(non_canonical_exons))

#re-open file and overwrite previous data. now will have exon IDs in correct order
file = open("crit_exon.txt","w")
for i in range(len(exon_locations)):
  file.write("{} = {}\n".format(result[i],exon_locations[i]))
file.close()

#write non_canonical_exon IDs into a file
file = open("non_canonical_exons.txt","w")
for i in range(len(non_canonical_exons)):
  file.write("{}\n".format(non_canonical_exons[i]))
file.close()

file = open("final.txt","w")
descount = 0
#write file including target genes, transcript IDs, exon IDs, phases, and rank
#looking for genes with single critical exon, and with multiple
for i in exon_locations:
    try:
        if(type(i) == list):
            file.write("{} = {},{},{},{},{},{}\n".format(gene_names[i[0]],trans_ids_exon_first[i[0]],i[0],i[1],exon_dict[i[0]][1],exon_dict[i[1]][2],rank[i[0]]))
            descount+=1
        else:
            #format gene name = transcript ID, crit exon ID, start phase, end phase, rank
            file.write("{} = {},{},{},{},{}\n".format(gene_names[i],trans_ids_exon_first[i],i,exon_dict[i][1],exon_dict[i][2],rank[i]))
            descount+=1
    except KeyError:
        continue

file.close()

print("Number of exons ready for analysis: ", descount)

#list to store start and end phases for each gene
gene_phases = {}
no_coding = [] # FORMAT = [ exonid, geneid, gene name ]

#CHECK if exon ID exists or not in pcg canonical file
def check(num):
    if (num not in exon_dict.keys()):
        no_coding.append(num)
        return False
    else:
        return True

#FORMAT
# exon_dict - [ transid, start phase , end phase ]
exon_mismatch=[]
#for loop looking for the ID of the exon that comes before the critical exon (rank-1)
#for genes with multiple critical exon, the one with lower rank is used as starting point
#rank-1 exon IDs uploaded to BioMart to obtain sequences
#check the multi-exon designs to make sure they are from the same gene - if they are not (which is the case for Stfa2 and Stfa2l, and it's a problem with the source data), ignore that entry
for i in exon_locations:
    if (type(i) == list):
        tupkey = tuple(i)
        if(not check(i[0]) or (not check(i[1]))): continue
        gene_ids[tupkey] = gene_ids[tupkey[0]]
        if any( gene_ids[x] != gene_ids[tupkey] for x in i ):
          exon_mismatch.append(tupkey)
        else:
          gene_phases[tupkey] = [exon_dict[i[0]][1],exon_dict[i[1]][-1]]
    else:
        tupkey=(i,)
        if (not check(i)): continue
        gene_ids[tupkey] = gene_ids[tupkey[0]]
        gene_phases[tupkey] = [exon_dict[i][1],exon_dict[i][-1]]

print("Number of designs with exons from different genes: ",len(exon_mismatch))
print("Number of designs ready for analysis: ", len(gene_phases.keys()))

exon_mismatch_file = open('exon_mismatch.txt','w')
for i in exon_mismatch:
    exon_mismatch_file.write("{}\n".format(i))
exon_mismatch_file.close()


#lists for different categories
stop_codon_en2 = [] #phase start 2
read_through = [] #phase 0-1 OR 1-2
frameshift = [] #phase 0-0 OR 0-2 OR 1-0 OR 1-1
truncated_protein = [] #phase start 2, ends in  -1 (these are a subset of what is in stop_codon_en2)
start_phase_negative = [] #phase start -
end_phase_negative =[] #end with - (not including the genes in stop_codon_en2 and truncated_protein)
unaccounted = []

for i in gene_phases.keys():
    if (gene_phases[i][0] == 2):
        stop_codon_en2.append(i)
    elif (gene_phases[i][0] < 0):
        start_phase_negative.append(i)
    elif (gene_phases[i][1] < 0):
        end_phase_negative.append(i)
    elif ((gene_phases[i][0] == 0 and gene_phases[i][1] == 1)):
        read_through.append(i)
    elif ((gene_phases[i][0] == 1 and gene_phases[i][1] == 2)):
        read_through.append(i)
    elif ((gene_phases[i][0] == 0 and gene_phases[i][1] == 2)):
        frameshift.append(i)
    elif ((gene_phases[i][0] == 1 and gene_phases[i][1] == 1)):
        frameshift.append(i)
    elif ((gene_phases[i][0] == 0 and gene_phases[i][1] == 0)):
        frameshift.append(i)
    elif ((gene_phases[i][0] == 1 and gene_phases[i][1] == 0)):
        frameshift.append(i)
    if ((gene_phases[i][0] == 2 and gene_phases[i][1] == -1)):
        truncated_protein.append(i)

print("Stop_codon_en2:",len(stop_codon_en2))
print("Readthrough:",len(read_through))
print("Frameshift:",len(frameshift))
print("Unaccounted:",len(unaccounted))

#dictionaries of gene ids
duplicate_genes = {}
stop_ids = {}
fs_ids = {}
rt_ids = {}
ns_ids = {}
npend_ids = {}
trunc_ids = {}

stop_codon_file = open('stop_codon.txt','w')
for i in stop_codon_en2:
    stop_codon_file.write("{} = {},{}\n".format(i,gene_ids[i],gene_names[i[0]]))
    stop_ids[gene_ids[i]] = 1;
stop_codon_file.close()

#this is an extra collection of stop_codon genes which end with a negative phase, and thus should not be checked for duplicates - they have already been checked with the rest of the stop_codon genes
truncated_protein_file = open('truncated_protein.txt','w')
for i in truncated_protein:
    truncated_protein_file.write("{} = {},{}\n".format(i,gene_ids[i],gene_names[i[0]]))
    trunc_ids[gene_ids[i]] = 1;
truncated_protein_file.close()

read_through_file = open('read_through.txt','w')
for i in read_through:
    read_through_file.write("{} = {},{}\n".format(i,gene_ids[i],gene_names[i[0]]))
    if ((gene_ids[i] in stop_ids) and (gene_ids[i] not in rt_ids)):
      if (gene_ids[i] in duplicate_genes):
        duplicate_genes[gene_ids[i]] = duplicate_genes[gene_ids[i]]+1
      else:
        duplicate_genes[gene_ids[i]] = 2
    rt_ids[gene_ids[i]] = 1;
read_through_file.close()

frameshift_file = open('frameshift.txt','w')
for i in frameshift:
    frameshift_file.write("{} = {},{}\n".format(i,gene_ids[i],gene_names[i[0]]))
    if (((gene_ids[i] in stop_ids) or (gene_ids[i] in rt_ids)) and (gene_ids[i] not in fs_ids) ):
      if (gene_ids[i] in duplicate_genes):
        duplicate_genes[gene_ids[i]] = duplicate_genes[gene_ids[i]]+1
      else:
        duplicate_genes[gene_ids[i]] = 2
    fs_ids[gene_ids[i]] = 1;
frameshift_file.close()

start_phase_negative_file = open('start_phase_negative.txt','w')
for i in start_phase_negative:
    start_phase_negative_file.write("{} = {},{}\n".format(i,gene_ids[i],gene_names[i[0]]))
    if (((gene_ids[i] in stop_ids) or (gene_ids[i] in rt_ids) or (gene_ids[i] in fs_ids)) and (gene_ids[i] not in ns_ids)):
      if (gene_ids[i] in duplicate_genes):
        duplicate_genes[gene_ids[i]] = duplicate_genes[gene_ids[i]]+1
      else:
        duplicate_genes[gene_ids[i]] = 2
    ns_ids[gene_ids[i]] = 1;
start_phase_negative_file.close()

end_phase_negative_file = open('end_phase_negative.txt','w')
for i in end_phase_negative:
    end_phase_negative_file.write("{} = {},{}\n".format(i,gene_ids[i],gene_names[i[0]]))
    if (((gene_ids[i] in stop_ids) or (gene_ids[i] in rt_ids) or (gene_ids[i] in fs_ids) or (gene_ids[i] in ns_ids)) and (gene_ids[i] not in npend_ids)):
      if (gene_ids[i] in duplicate_genes):
        duplicate_genes[gene_ids[i]] = duplicate_genes[gene_ids[i]] + 1
      else:
        duplicate_genes[gene_ids[i]] = 2
    npend_ids[gene_ids[i]] = 1;
end_phase_negative_file.close()

print("Number of duplicate genes: ", len(duplicate_genes))

#go through the negative phase categories and count the duplicates
dup_count = 0
for i in npend_ids:
        if ( i in duplicate_genes ):
            dup_count+=1
print("End phase negative count (excluding duplicates): ", (len(npend_ids) - dup_count))

dup_count = 0
for i in ns_ids:
        if ( i in duplicate_genes ):
            dup_count+=1
print("Start phase negative count (excluding duplicates): ", (len(ns_ids) - dup_count))

dup_count = 0
for i in trunc_ids:
        if ( i in duplicate_genes ):
            dup_count+=1
print("Truncated count (excluding duplicates): ", (len(trunc_ids) - dup_count))


#output duplicate gene list
dup_gene_file = open('duplicate_designs.txt','w')
for i in duplicate_genes:
    dup_gene_file.write("{}\t{}\n".format(i,duplicate_genes[i]))
dup_gene_file.close()

#Parsing .rpt file as csv; file from https://www.informatics.jax.org/downloads/reports/index.html
filename = "MRK_ENSEMBL.rpt"

gene_model = {}

with open(filename, "r") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        # extract the data from each row
        gene_model[row[5]] = (row[1],row[0].split(":")[1])

file = open('mgi_match.txt', 'w')
missing = []
for i in gene_phases.keys():
    try:
        file.write("{} = {},{}\n".format(gene_model[gene_ids[i]][0],gene_ids[i],gene_model[gene_ids[i]][1]))
    except KeyError:
        missing.append(i)

file.close()


gene_dict = {}
gene_dict_del = {}

#parse procedure completeness and phenotype hits file. create two dictionaries:
#gene_dict for alleles where critical exon still present: tm1a, tm1e, tm1e.1
#gene_dict_del for alleles where exon is gone: tm1, tm1.1, tm1b
#file from http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/
with open('procedureCompletenessAndPhenotypeHits.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        geneid = row['Gene Accession Id'].split(":")[1].strip("\"")
        gene_symbol = row['Gene Symbol'].strip("\"")
        allele_symbol = row['Allele Symbol'].strip("\"")
        mp_term_id = row['MP Accession Id: Significant'].strip("\"")
        phenotype = row['MP Name: Significant'].strip("\"")
        num_p = row['MP Count: Significant'].strip("\"")
        if (re.search('<tm\d[ae](.1)?\(', allele_symbol)):
            if (geneid in gene_dict):
                if (allele_symbol not in gene_dict[geneid][0]):
                    new_a = gene_dict[geneid][0] + "|" + allele_symbol
                    gene_dict[0] = new_a
                if (mp_term_id not in gene_dict[geneid][1]):
                    new_p = gene_dict[geneid][2] + "|" + phenotype
                    new_t = gene_dict[geneid][1] + "|" + mp_term_id
                    new_n = int(gene_dict[geneid][3]) + int(num_p)
                    gene_dict[geneid] = (gene_dict[geneid][0], new_t, new_p, new_n)
            else:
                gene_dict[geneid] = (allele_symbol, mp_term_id, phenotype, num_p)
        elif (re.search('<tm\db?(.1)?\(', allele_symbol)):
            if (geneid in gene_dict_del):
                if (allele_symbol not in gene_dict_del[geneid][0]):
                    new_a = gene_dict_del[geneid][0] + "|" + allele_symbol
                    gene_dict_del[0] = new_a
                if (mp_term_id not in gene_dict_del[geneid][1]):
                    new_p = gene_dict_del[geneid][2] + "|" + phenotype
                    new_t = gene_dict_del[geneid][1] + "|" + mp_term_id
                    new_n = int(gene_dict_del[geneid][3]) + int(num_p)
                    gene_dict_del[geneid] = (gene_dict_del[geneid][0], new_t, new_p, new_n)
            else:
                gene_dict_del[geneid] = (allele_symbol, mp_term_id, phenotype, num_p)


#example
print(list(gene_dict.values())[3])
print(list(gene_dict_del.values())[3])

#go through 3 main categories, using the gene IDs because the phenotypes are counted by gene, not by design

#stop codon
file = open('stop_codon_phen.txt','w')
#list to store gene IDs that are within the IMPC targeted mutation file, but not tested for phenotypes
stop_codon_untested = []
#to count how many genes have 0 phenotypes (tested, but 0)
zero_count_ko = 0
zero_count_del = 0
#sum of phenotypes from knockout-first and deletion alleles
ko_sum = 0
del_sum = 0
#sum of lines assessed (by category - if there was a tm1a and a tm1e created for the same gene they would be included as 1 count)
ko_lines = 0
del_lines = 0
all_lines = 0
#sum of preweaning lethalities
lethality_ko_count = 0
lethality_del_count = 0
#sum of duplicates
dup_count = 0
for i in stop_ids:
    #exclude duplicates
    if ( i in duplicate_genes ):
        dup_count+=1
    #write into a file gene name and no. of phenotypes
    elif ((gene_model[i][1] in gene_dict) or (gene_model[i][1] in gene_dict_del)):
        all_lines +=1
        if (gene_model[i][1] in gene_dict):
            file.write("{},{} = {},{}\n".format(i,gene_dict[gene_model[i][1]][0],gene_dict[gene_model[i][1]][3],gene_dict[gene_model[i][1]][2]))
            ko_sum += int(gene_dict[gene_model[i][1]][3])
            if ("lethality" in gene_dict[gene_model[i][1]][2]):
                lethality_ko_count += 1
            ko_lines += 1
            if (int(gene_dict[gene_model[i][1]][3]) == 0):
              zero_count_ko += 1
        if (gene_model[i][1] in gene_dict_del):
            file.write("{},{} = {},{}\n".format(i,gene_dict_del[gene_model[i][1]][0],gene_dict_del[gene_model[i][1]][3],gene_dict_del[gene_model[i][1]][2]))
            del_sum += int(gene_dict_del[gene_model[i][1]][3])
            if ("lethality" in gene_dict_del[gene_model[i][1]][2]):
                lethality_del_count += 1
            del_lines += 1
            if (int(gene_dict_del[gene_model[i][1]][3]) == 0):
              zero_count_del += 1
    #write in the gene names that are not tested into list
    else:
        stop_codon_untested.append(i)
file.close()

print("Stop Codon design count:", len(stop_codon_en2))
print("Stop Codon gene count:", len(stop_ids))
print("Stop Codon gene count (excluding duplicates):", (len(stop_ids) - dup_count))
print("duplicates (excluded):",dup_count)
print("Stop_Codon tested",all_lines)
print("Stop_Codon Untested",len(stop_codon_untested))
print("Knockout-first (tm1a, tm1e, tm1e.1):")
print("Lines with no phenotypes: ", zero_count_ko)
print("Sum of phenotypes: ", ko_sum)
print("Lines with lethality: ",lethality_ko_count)
print("Lines tested: ",ko_lines)
print("Average phenotype count/tested line: ", ko_sum/ko_lines)
print("Reporter-tagged deletion alleles (tm1b, tm1, tm1.1):")
print("Lines with no phenotypes: ", zero_count_del)
print("Sum of phenotypes: ", del_sum)
print("Lines with lethality: ",lethality_del_count)
print("Lines tested: ",del_lines)
print("Average phenotype count/tested line: ", del_sum/del_lines)

#Readthrough
file = open('read_through_phen.txt','w')
read_through_untested = []
zero_count_ko = 0
zero_count_del = 0
ko_sum = 0
del_sum = 0
ko_lines = 0
del_lines = 0
all_lines = 0
lethality_ko_count = 0
lethality_del_count = 0
dup_count = 0
for i in rt_ids:
    #exclude duplicates
    if ( i in duplicate_genes ):
        dup_count+=1
    #write into a file gene name and no. of phenotypes
    elif ((gene_model[i][1] in gene_dict) or (gene_model[i][1] in gene_dict_del)):
        all_lines +=1
        if (gene_model[i][1] in gene_dict):
            file.write("{},{} = {},{}\n".format(i,gene_dict[gene_model[i][1]][0],gene_dict[gene_model[i][1]][3],gene_dict[gene_model[i][1]][2]))
            ko_sum += int(gene_dict[gene_model[i][1]][3])
            if ("lethality" in gene_dict[gene_model[i][1]][2]):
                lethality_ko_count += 1
            ko_lines += 1
            if (int(gene_dict[gene_model[i][1]][3]) == 0):
              zero_count_ko += 1
        if (gene_model[i][1] in gene_dict_del):
            file.write("{},{} = {},{}\n".format(i,gene_dict_del[gene_model[i][1]][0],gene_dict_del[gene_model[i][1]][3],gene_dict_del[gene_model[i][1]][2]))
            del_sum += int(gene_dict_del[gene_model[i][1]][3])
            if ("lethality" in gene_dict_del[gene_model[i][1]][2]):
                lethality_del_count += 1
            del_lines += 1
            if (int(gene_dict_del[gene_model[i][1]][3]) == 0):
              zero_count_del += 1
    #write in the gene names that are not tested into list
    else:
        read_through_untested.append(i)
file.close()
print("Readthrough design count:", len(read_through))
print("Readthrough gene count:", len(rt_ids))
print("Readthrough gene count (excluding duplicates):", (len(rt_ids) - dup_count))
print("duplicates (excluded):",dup_count)
print("Readthrough tested:",all_lines)
print("Readthrough Untested:",len(read_through_untested))
print("Knockout-first (tm1a, tm1e, tm1e.1):")
print("Lines with no phenotypes: ", zero_count_ko)
print("Sum of phenotypes: ", ko_sum)
print("Lines with lethality: ",lethality_ko_count)
print("Lines tested: ",ko_lines)
print("Average phenotype count/tested line: ", ko_sum/ko_lines)
print("Reporter-tagged deletion alleles (tm1b, tm1, tm1.1):")
print("Lines with no phenotypes: ", zero_count_del)
print("Sum of phenotypes: ", del_sum)
print("Lines with lethality: ",lethality_del_count)
print("Lines tested: ",del_lines)
print("Average phenotype count/tested line: ", del_sum/del_lines)


#Frameshift
file = open('frameshift_phen.txt','w')
frameshift_untested = []
zero_count_ko = 0
zero_count_del = 0
ko_sum = 0
del_sum = 0
ko_lines = 0
del_lines = 0
all_lines = 0
lethality_ko_count = 0
lethality_del_count = 0
dup_count = 0
for i in fs_ids:
    #exclude duplicates
    if ( i in duplicate_genes ):
        dup_count+=1
    #write into a file gene name and no. of phenotypes
    elif ((gene_model[i][1] in gene_dict) or (gene_model[i][1] in gene_dict_del)):
        all_lines +=1
        if (gene_model[i][1] in gene_dict):
            file.write("{},{} = {},{}\n".format(i,gene_dict[gene_model[i][1]][0],gene_dict[gene_model[i][1]][3],gene_dict[gene_model[i][1]][2]))
            ko_sum += int(gene_dict[gene_model[i][1]][3])
            if ("lethality" in gene_dict[gene_model[i][1]][2]):
                lethality_ko_count += 1
            ko_lines += 1
            if (int(gene_dict[gene_model[i][1]][3]) == 0):
              zero_count_ko += 1
        if (gene_model[i][1] in gene_dict_del):
            file.write("{},{} = {},{}\n".format(i,gene_dict_del[gene_model[i][1]][0],gene_dict_del[gene_model[i][1]][3],gene_dict_del[gene_model[i][1]][2]))
            del_sum += int(gene_dict_del[gene_model[i][1]][3])
            if ("lethality" in gene_dict_del[gene_model[i][1]][2]):
                lethality_del_count += 1
            del_lines += 1
            if (int(gene_dict_del[gene_model[i][1]][3]) == 0):
              zero_count_del += 1
    #write in the gene names that are not tested into list
    else:
        frameshift_untested.append(i)

file.close()
print("Frameshift design count:", len(frameshift))
print("Frameshift gene count:", len(fs_ids))
print("Frameshift gene count (excluding duplicates):", (len(fs_ids) - dup_count))
print("duplicates (excluded):",dup_count)
print("Frameshift tested:",all_lines)
print("Frameshift Untested:",len(frameshift_untested))
print("Knockout-first (tm1a, tm1e, tm1e.1):")
print("Lines with no phenotypes: ", zero_count_ko)
print("Sum of phenotypes: ", ko_sum)
print("Lines with lethality: ",lethality_ko_count)
print("Lines tested: ",ko_lines)
print("Average phenotype count/tested line: ", ko_sum/ko_lines)
print("Reporter-tagged deletion alleles (tm1b, tm1, tm1.1):")
print("Lines with no phenotypes: ", zero_count_del)
print("Sum of phenotypes: ", del_sum)
print("Lines with lethality: ",lethality_del_count)
print("Lines tested: ",del_lines)
print("Average phenotype count/tested line: ", del_sum/del_lines)


