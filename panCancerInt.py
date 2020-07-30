#!/usr/bin/env python3
#./pythonScript.py --i testCosmicBreastReccurentMutationsAnnotated.bed  --o test.tmp --r ../../RNA-Seq/gene_expression_tophatStarFPKM_v2.tsv --t ../../breast/breastPatientIds.txt
#Wrote a file called IdConversion.tsv with the biomart gene symbol conversions and used R to edit the RNA-Seq table 
#RNASeq<- read.table("../../RNA-Seq/gene_expression_tophatStarFPKM_v2.tsv")
#ids<-read.delim("IdConversion.tsv",header=F)
#write.csv(cbind(ids,RNASeq), "RNA-Seq.csv")
# sed 1d RNA-Seq.csv -i

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Take the file CosmicBreastReccurentMutationsAnnotated.bed and try to find mutations linked to differential expression.")
parser.add_argument("--input",help="Input file (*ReccurentMutationsAnnotated.bed).",required=True)
parser.add_argument("--output",help="Output file",required=True)
parser.add_argument("--rnaseq",help="RNA-Seq quant file location",required=True)
parser.add_argument("--tissueSamples",help="Patients affected by the target cancer type",required=True)

args=parser.parse_args()

f= open(args.output,'w')

#Create a dictionary
d={}

#Read list of patients with cancer type
with open(args.tissueSamples) as t:
    patients = t.readlines()
patients = [x.strip() for x in patients]
#print(patients)

#Read the RNA-Seq file into a pandas df
rnaDF=pd.read_csv(args.rnaseq)

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

#Get a list of the patients which have rna-seq data
allPatients=list(rnaDF.columns)
#print(allPatients)
patients=intersection(allPatients, patients)
patients.insert(0,"feature")
patients.insert(0,"geneName")
patients.insert(0,"geneid")
patients.insert(0,"1")
print(patients)
#subset to just those with affected tissue type
rnaDF=rnaDF[patients]
#print(rnaDF.iloc[1:10,"geneName"])

#rnaDF2=rnaDF.iloc[["geneName","feature","01370d42-f75c-4532-9b9c-24ff7302b033"]]
#print(rnaDF2.head())
#geneid geneName             feature  01370d42-f75c-4532-9b9c-24ff7302b033
#loop through each mutation
for line in open(args.input):
    #Get the donorID hashes and store them into files
    files=line.split("\t")[1]
    #Strip newline from gene name to be used in pandas slice below
    gene=line.split("\t")[6].strip()
    mutation=line.split("\t")[0]

    donors=files.split(",")
    print(donors)
    donors=intersection(donors, patients)
    #This line works, and will select the correct cols
    #print(rnaDF.loc[rnaDF["geneName"] == gene,["geneid","01370d42-f75c-4532-9b9c-24ff7302b033"]])
    #print(rnaDF.loc[[rnaDF['geneName'] == gene],["geneid","01370d42-f75c-4532-9b9c-24ff7302b033"]])
    #print(rnaDF.loc[[rnaDF['geneName'] == gene],[donors]])
#Loop through each file and get the quantification of the gene 
