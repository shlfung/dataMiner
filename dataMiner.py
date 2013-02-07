# By Stephen Fung
# The script reads the OFP result and Automated ID result based on a given rank and a name. It outputs the the data into a CSV file


import csv
import re
import sys
import glob
import os
import math
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter
from collections import namedtuple
from datetime import datetime


# The metaReader function reads the metaData File and return multiple dictionaries:
# CSV File --> Dictionaries
def metaReader(metaData):

    #Begin Reading CSV File
    metaDataFile = open(metaData, "rb")
    metaDataReader = csv.reader(metaDataFile, delimiter=",")
    metaDataHeader = metaDataReader.next() #Skip the Header

    SFF_MID_2_SampleID = {} # SFF_MID_2_Sample for finding out Sample ID using SFF.MID as keys
    SampleID_2_Collector = {}   # SampleID_2_Collector for finding out Collector using Sample ID as keys
    SampleID_2_Year = {}    # SampleID_2_Year for finding out Years using Sample ID as keys
    SampleID_2_Week = {}    # SampleID_2_Week for finding out Weeks using Sample ID as keys
    SampleID_2_City = {}    # SampleID_2_City for finding out Cities using Sample ID as keys
    SampleID_2_Province = {}    # SampleID_2_Province for finding out Provinces using Sample ID as keys
    SampleID_2_Dates = {}   # SampleID_2_Dates for actual calendar dates using Sample ID as keys
    SampleID_2_MaxTemp = {} # SampleID_2_MaxTemp for maximum daily temperature using Sample ID as keys
    SampleID_2_MinTemp = {} # SampleID_2_MinTemp for minimum daily temperature using Sample ID as keys
    SampleID_2_MeanTemp = {}  # SampleID_2_MeanTemp for average daily temperature using Sample ID as keys
    SampleID_2_Percipitation = {}  # SampleID_2_Percipitation for daily total percipitation using Sample ID as keys
    SampleID_2_WindDir = {}   # SampleID_2_WindDir for the wind direction using Sample ID as keys
    SampleID_2_MaxWindSpd = {}  # SampleID_2_MaxWindSpd for daily maximum wind speed using Sample ID as keys

    
    
    #Iterate over the CSV file, read it row by row, and store the information to the appropriate dictionaries
    for row in metaDataReader:
        SFF_MID = str(row[3]) + "." + 'MID' + str(row[6])
        SFF_MID_2_SampleID[SFF_MID] = row[7] #Generate the SFF.MID using row[3] and row[6]
        SampleID_2_Collector[row[7]] = row[8]
        SampleID_2_Year[row[7]] = row[9]
        SampleID_2_Week[row[7]] = row[10]
        SampleID_2_City[row[7]] = row[13]
        SampleID_2_Province[row[7]] = row[14]
        try:
            SampleID_2_Dates[row[7]] = str(row[33]) + '-' + str(row[9])
        except IndexError:
            pass
        try:
            SampleID_2_MaxTemp[row[7]] = row[34]
        except IndexError:
            pass
        try:
            SampleID_2_MinTemp[row[7]] = row[35]
        except IndexError:
            pass
        try:
            SampleID_2_MeanTemp[row[7]] = row[36]
        except IndexError:
            pass
        try:
            SampleID_2_Percipitation[row[7]] = row[41]
        except IndexError:
            pass
        try:
            SampleID_2_WindDir[row[7]] = row[43]
        except IndexError:
            pass
        try:
            SampleID_2_MaxWindSpd[row[7]] = row[44]
        except IndexError:
            pass
        
    return SFF_MID_2_SampleID, SampleID_2_Collector, SampleID_2_Year, SampleID_2_Week, SampleID_2_City, SampleID_2_Province, SampleID_2_Dates, SampleID_2_MaxTemp, SampleID_2_MinTemp, SampleID_2_MeanTemp, SampleID_2_Percipitation, SampleID_2_WindDir, SampleID_2_MaxWindSpd

#The summary_log_reader function reads alll the summary log files row by row,
# When it finds a row with the matching rank and name, it stores the name of the SFF.MID number to a list
# Summary Log Files, String, String --> List
def summary_log_reader(summarylogs, name, rank):
    
    summarylogs = summarylogs + '*.summarylog'
    listOfSummaryLogs = glob.glob(summarylogs) # Using glob to read multiple files in a directory
    allSummaryLogs = []

    # Read each file of the summary log
    for aFile in listOfSummaryLogs:
        fileName = str(aFile)
        print 'Reading', fileName
        sep = '.'
        MIDNumber = fileName.split(sep)[2]
        #MID = re.sub("[^0-9]", "", MID)
        MIDNumber = str(MIDNumber).strip()   #Extract the SFF.MID number from the file name
        aSummaryLogFile = open(aFile, "rb")
        aReader = csv.reader(aSummaryLogFile, delimiter = "\t")
        # Check if each row in a file has a matching name and rank with the input
        for row in aReader:
            if rank == 'kingdom':
                try:
                    if row[13] == name:
                        aName = str(row[0]).strip() + '.' + MIDNumber
                        allSummaryLogs.append(aName)
                        print 'Found', row[13]
                        print 'Found', len(allSummaryLogs), 'sequences'
                    else:
                        pass
                except IndexError:
                    pass
            if rank == 'phylum':
                try:
                    if row[14] == name:
                        aName = str(row[0]).strip() + '.' + MIDNumber
                        allSummaryLogs.append(aName)
                        print 'Found', row[14]
                        print 'Found', len(allSummaryLogs), 'sequences'
                    else:
                        pass
                except IndexError:
                    pass
            if rank == 'class':
                try:
                    if row[15] == name:
                        aName = str(row[0]).strip() + '.' + MIDNumber
                        allSummaryLogs.append(aName)
                        print 'Found', row[15]
                        print 'Found', len(allSummaryLogs), 'sequences'
                    else:
                        pass
                except IndexError:
                    pass
            if rank == 'order':
                try:
                    if row[16] == name:
                        aName = str(row[0]).strip() + '.' + MIDNumber
                        allSummaryLogs.append(aName)
                        print 'Found', row[16]
                        print 'Found', len(allSummaryLogs), 'sequences'
                    else:
                        pass
                except IndexError:
                    pass
            if rank == 'family':
                try:
                    if row[17] == name:
                        aName = str(row[0]).strip() + '.' + MIDNumber
                        allSummaryLogs.append(aName)
                        print 'Found', row[17]
                        print 'Found', len(allSummaryLogs), 'sequences'
                    else:
                        pass
                except IndexError:
                    pass
            if rank == 'genus':
                try:
                    if row[18] == name:
                        aName = str(row[0]).strip() + '.' + MIDNumber
                        allSummaryLogs.append(aName)
                        print 'Found', row[18]
                        print 'Found', len(allSummaryLogs), 'sequences'
                    else:
                        pass
                except IndexError:
                    pass
            if rank == 'species':
                try:
                    sep = ' '
                    genus = (row[19].split(sep)[0]).strip()
                    species = (row[19].split(sep)[1]).strip()
                    theNameToCheck = genus + ' ' + species
                    if name == theNameToCheck:
                        aName = str(row[0]).strip() + '.' + MIDNumber
                        allSummaryLogs.append(aName)
                        print 'Found', theNameToCheck
                        print 'Found', len(allSummaryLogs), 'sequences'
                    else:
                        pass
                except IndexError:
                    pass
    return allSummaryLogs 
        
# The oligoFishingReader function reads the output from the Oligo Fishing Pipeline (OFP) and stores the Sample ID number of each sequence to a list
# Text in fasta format --> List
def oligoFishingReader(ofrFile):
    
    seqRecord = SeqIO.parse(ofrFile, "fasta") #Using Biopython module to parse the fasta file

    listOfSampleIDs = []

    # Read each fasta sequence
    for i in seqRecord:
        descr = i.description
        aSequenceProperty = dict(x.split(":") for x in re.findall(r"\w+:\s\w+", descr))
        # Each fasta sequences has multiple dscriptions. This breaks down the description into multiple categories using regular expression.
        # The categories are broken up into dictionaries temporarily 
        if aSequenceProperty['Project'] == ' SAGES':
            listOfSampleIDs.append(aSequenceProperty['Sample'])
        else:
            continue

    listOfSampleIDs = [i.strip() for i in listOfSampleIDs] #Remove empty spaces from the Sample IDs
    
    return listOfSampleIDs

# This function removes keys of a dictionary that are not present in the lists
# Dict, List, List -> Dict
def removeRedundantMeta(aDict, aList, aList2):

    aSet = set(aList)
    aSet2 = set(aList2)
    
    for k in aDict.keys():
        if (k not in aSet) and (k not in aSet2):
            del aDict[k]

    return aDict
    
# This is the main function that calls the other functions
# It produces a CSV file that compare the oligo fishing and automated identification result.
# It takes the rank, name, the summary logs, the metaData file, the OFP output
# String, string, summary log files, metaData File, OFP fasta text file --> CSV file or DataFrame
def analyzer(rank, name, path_to_summary_logs, path_to_metaData, path_to_ofr_result, removeRedundant):

    startTime = datetime.now()
    print 'Now opening all the summary logs and store to memory'
    
    summaryLogs = summary_log_reader(path_to_summary_logs, name, rank)
    
    listOfOrfElements = oligoFishingReader(path_to_ofr_result)
    
    multipleDicts = metaReader(path_to_metaData) 
    metaSampleIDs = multipleDicts[0]
    metaCollectors = multipleDicts[1]
    metaYears = multipleDicts[2]
    metaWeeks = multipleDicts[3]
    metaCities = multipleDicts[4]
    metaProvinces = multipleDicts[5]
    metaDates = multipleDicts[6]
    metaMaxTemp = multipleDicts[7]
    metaMinTemp = multipleDicts[8]
    metaMeanTemp = multipleDicts[9]
    metaPercipitation = multipleDicts[10]
    metaWindDir = multipleDicts[11]
    metaMaxWindSpd = multipleDicts[12]

    # Here we convert the list of SFF.MID from the summary_log_reader function to a list of Sample IDs
    listOfAutoElements = []
    for i in summaryLogs:
        #print row
        try: 
            aSampleID = metaSampleIDs[i]
            listOfAutoElements.append(aSampleID)   
        except KeyError:
            pass

    # Remove the extra metaData information that are not present from the automated ID data or the oligo fishing data
    if removeRedundant == 'Yes':
        metaCollectors = removeRedundantMeta(metaCollectors, listOfAutoElements, listOfOrfElements)
        metaYears = removeRedundantMeta(metaYears, listOfAutoElements, listOfOrfElements)
        metaWeeks = removeRedundantMeta(metaWeeks, listOfAutoElements, listOfOrfElements)
        metaCities = removeRedundantMeta(metaCities, listOfAutoElements, listOfOrfElements)
        metaProvinces = removeRedundantMeta(metaProvinces, listOfAutoElements, listOfOrfElements)
        metaDates = removeRedundantMeta(metaDates, listOfAutoElements, listOfOrfElements)
        metaMaxTemp = removeRedundantMeta(metaMaxTemp, listOfAutoElements, listOfOrfElements)
        metaMinTemp = removeRedundantMeta(metaMinTemp, listOfAutoElements, listOfOrfElements)
        metaMeanTemp = removeRedundantMeta(metaMeanTemp, listOfAutoElements, listOfOrfElements)
        metaPercipitation = removeRedundantMeta(metaPercipitation, listOfAutoElements, listOfOrfElements)
        metaWindDir = removeRedundantMeta(metaWindDir, listOfAutoElements, listOfOrfElements)
        metaMaxWindSpd = removeRedundantMeta(metaMaxWindSpd, listOfAutoElements, listOfOrfElements)
    else:
        pass

    # Count the appearance of each element in the list of samples IDs and output the results into dictionaries
    numberOfOrfElements = Counter(listOfOrfElements)
    numberOfAutoElements = Counter(listOfAutoElements)
    
    # Converting all the data into a dataframe using pandas
    print 'Generating DataFrame and CSV output...'
    temp = {'Provinces': pd.Series(metaProvinces), 'Cities': pd.Series(metaCities), 'Collectors': pd.Series(metaCollectors), 'Years': pd.Series(metaYears), 'Weeks': pd.Series(metaWeeks), 'Dates': pd.Series(metaDates), 'Max Temperature': pd.Series(metaMaxTemp), 'Min Temperature': pd.Series(metaMinTemp), 'Mean Temperature': pd.Series(metaMeanTemp), 'Total Percipitation': pd.Series(metaPercipitation), 'Wind Direction': pd.Series(metaWindDir), 'Max Wind Speed': pd.Series(metaMaxWindSpd), 'Auto ID': pd.Series(numberOfAutoElements), 'OFP Result': pd.Series(numberOfOrfElements)}
    df = pd.DataFrame(temp, columns=['Provinces', 'Cities', 'Collectors', 'Years', 'Weeks', 'Dates', 'Max Temperature', 'Min Temperature', 'Mean Temperature', 'Total Percipitation', 'Wind Direction', 'Max Wind Speed', 'Auto ID', 'OFP Result'])
    df = df.fillna(value=0)
    df['Log(Auto ID + 1)'] = df['Auto ID'].apply(np.log1p)
    df['Log(OFP Result + 1)'] = df['OFP Result'].apply(np.log1p)
    fileName = rank + '_' + name + '_' + 'Auto vs OFP.csv'
    df.to_csv(fileName)
    print 'Took', (datetime.now()-startTime)
    return df
        
    

if __name__=='__main__':
    #temp = metaReader('./input/454_sample_summary_with_climate.csv')
    #temp2 = summary_log_reader('./rs75 Summary Log/', 'Staphylococcus aureus', 'species')
    #temp3 = oligoFishingReader('./input/Streptococcus_spp.Streptococcus_spp.txt.withmeta.fna.fasta')
##    item1 = analyzer('species', 'Staphylococcus aureus', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Staphylococcus_aureus.Staphylococcus_aureus.txt.withmeta.fna.fasta', 'No')
##    item2 = analyzer('species', 'Bacillus anthracis', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Bacillus_anthracis_cand1.Bacillus_anthracis_cand1.txt.withmeta.fna.fasta', 'No')
##    item3 = analyzer('species', 'Bacillus herbersteinensis', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Bacillus_anthracis_cand1.Bacillus_anthracis_cand1.txt.withmeta.fna.fasta', 'No')
##    item4 = analyzer('species', 'Staphylococcus epidermidis', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Staphylococcus_epidermidis.Staphylococcus_epidermidis.txt.withmeta.fna.fasta', 'No')
##    item5 = analyzer('species', 'Clostridium botulinum', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Clostridium_botulinum.Clostridium_botulinum.txt.withmeta.fna.fasta', 'No')
##    item6 = analyzer('species', 'Clostridium perfringens', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Clostridium_perfringens.Clostridium_perfringens.txt.withmeta.fna.fasta', 'No')
##    item7 = analyzer('phylum', 'Firmicutes', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Firmicutes_oligo1.Firmicutes.txt.withmeta.fna.fasta', 'No')
##    item8 = analyzer('genus', 'Mycobacterium', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Mycobacterium_spp.Mycobacterium_spp.txt.withmeta.fna.fasta', 'No')
##    item9 = analyzer('genus', 'Streptococcus', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Streptococcus_spp.Streptococcus_spp.txt.withmeta.fna.fasta', 'No')
##    item10 = analyzer('genus', 'Bacillus', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Bacillus_spp.Bacillus_spp.txt.withmeta.fna.fasta', 'No')
##    item11 = analyzer('genus', 'Staphylococcus', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Staphylococcus_spp.Staphylococcus_spp.txt.withmeta.fna.fasta', 'No')
##    item12 = analyzer('genus', 'Rickettsia', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Rickettsia_spp.Rickettsia_spp.txt.withmeta.fna.fasta', 'No')
    item13 = analyzer('genus', 'Clostridium', './rs75 Summary Log/', './input/454_sample_summary_with_climate.csv', './input/Clostridium_spp.Clostridium_spp.txt.withmeta.fna.fasta', 'No')




