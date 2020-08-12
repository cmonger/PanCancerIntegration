#!/usr/bin/env python3
#./pythonScript.py --i testCosmicBreastReccurentMutationsAnnotated.bed  --o test.tmp --r ../../RNA-Seq/gene_expression_tophatStarFPKM_v2.tsv --t ../../breast/breastPatientIds.txt
#Wrote a file called IdConversion.tsv with the biomart gene symbol conversions and used R to edit the RNA-Seq table 
#RNASeq<- read.table("../../RNA-Seq/gene_expression_tophatStarFPKM_v2.tsv")
#ids<-read.delim("IdConversion.tsv",header=F)
#write.csv(cbind(ids,RNASeq), "RNA-Seq.csv")
# sed 1d RNA-Seq.csv -i
#./panCancerInt.py --i CosmicBreastReccurentMutationsAnnotated.bed  --o test.tmp --r RNA-Seq.csv --t ../../breast/breastPatientIds.txt --s ../../pcawg_sample_sheet.tsv

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Take the file CosmicBreastReccurentMutationsAnnotated.bed and try to find mutations linked to differential expression.")
parser.add_argument("--input",help="Input file (*ReccurentMutationsAnnotated.bed).",required=True)
parser.add_argument("--output",help="Output file",required=True)
parser.add_argument("--rnaseq",help="RNA-Seq quant file location",required=True)
parser.add_argument("--tissueSamples",help="Patients affected by the target cancer type",required=True)
parser.add_argument("--sampleSheet",help="PanCancer sample sheet location",required=True)
args=parser.parse_args()

f= open(args.output,'w')


samples = pd.read_csv(args.sampleSheet,sep="\t")
#Use the aliquot id to get the donor id

#Function to convert WGS sample IDs to RNA-Seq ids
def convert(identifier):
    donorId=samples.loc[samples["aliquot_id"] == identifier,"donor_unique_id"].values[0]
    donorId2=samples.loc[(samples["library_strategy"] == "RNA-Seq")& (samples["donor_unique_id"] == donorId),"aliquot_id"]
    if donorId2.empty:
        return ""
    donorId2=samples.loc[(samples["library_strategy"] == "RNA-Seq")& (samples["donor_unique_id"] == donorId),"aliquot_id"].values[0]
    return donorId2


#Read list of patients with cancer type
with open(args.tissueSamples) as t:
    patients = t.readlines()
patients = [x.strip() for x in patients]
#print(patients)

#convert these to RNA-Seq aliquot ids and drop empty strings
patIds=[]
for pat in patients:
    patIds.append(convert(pat))
patIds = list(filter(None, patIds))
#print(patIds)

#Read the RNA-Seq file into a pandas df
rnaDF=pd.read_csv(args.rnaseq)

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def outersection(lst1,lst2):
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3

#Get a list of the patients which have rna-seq data
allPatients=list(rnaDF.columns)
#print(allPatients)
patients=intersection(allPatients, patIds)
patients.insert(0,"feature")
patients.insert(0,"geneName")
patients.insert(0,"geneid")
patients.insert(0,"1")
#print(patients)
#subset to just those with affected tissue type
rnaDF=rnaDF[patients]
#print (rnaDF)


#print(rnaDF.loc[1:10,"128ec04c-6514-4eb9-aa2e-79da6949b9a8"]) #Doesnt work
#print(rnaDF.iloc[1:10,"geneName"]) #Doesnt work
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
    #print(donors)
    RNAids= []
    for donorId in donors:
        RNAids.append(convert(donorId))
    #    print(convert(donorId))
    #print(RNAids)
    donors=intersection(RNAids, patIds)
    donors.insert(0,"geneName")
    #print(donors)
    #This line works, and will select the correct cols
    #print(rnaDF.loc[rnaDF["geneName"] == gene,["geneid","01370d42-f75c-4532-9b9c-24ff7302b033"]])
    #print(rnaDF.loc[[rnaDF['geneName'] == gene],["geneid","01370d42-f75c-4532-9b9c-24ff7302b033"]])
    #print(rnaDF.iloc[[rnaDF['geneName'] == gene],[donors]])
    
    if (len(donors) > 1):
        #Create a dataframe for patients RNA-Seq data with and without the mutation
        mutatedRNADF=rnaDF[rnaDF.columns[rnaDF.columns.isin(donors)]]
        nonMutatedIds=outersection(patIds,RNAids)
        nonMutatedRNADF=rnaDF[rnaDF.columns[rnaDF.columns.isin(nonMutatedIds)]]

        print(mutatedRNADF.loc[rnaDF["geneName"] == gene]) 
        print (nonMutatedRNADF.loc[rnaDF["geneName"] == gene])
#Loop through each file and get the quantification of the gene 


