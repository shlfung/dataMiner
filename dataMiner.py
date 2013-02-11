# By Stephen Fung
# The script reads the Oligo Fishing Pipeline (OFP) result and Automated Identification (Auto ID) result based on a given rank and a name.
# It outputs the Auto ID and OFP data and return the climate, geography...etc in a csv file that can be opened in MS Excel by the biologists
# The script requires scipy, numpy, biopython and pandas to be installed with python 2.7


import csv
import re
import sys
import glob
import os
import math
import collections
#from collections import Counter
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO



# The metaReader function reads the metaData File and return multiple dictionaries:
# The dictionaries uses Sample ID as keys for the relevant information
# CSV File --> Dictionaries
def meta_reader(meta_data):

    sff_mid_2_sample_id = {} # sff_mid_2_Sample for finding out Sample ID using SFF.MID as keys
    sample_id_2_collector = {}   # sample_id_2_collector for finding out Collector using Sample ID as keys
    sample_id_2_year = {}    # sample_id_2_year for finding out Years using Sample ID as keys
    sample_id_2_week = {}    # sample_id_2_week for finding out Weeks using Sample ID as keys
    sample_id_2_city = {}    # sample_id_2_city for finding out Cities using Sample ID as keys
    sample_id_2_province = {}    # sample_id_2_province for finding out Provinces using Sample ID as keys
    sample_id_2_dates = {}   # sample_id_2_dates for actual calendar dates using Sample ID as keys
    sample_id_2_max_temp = {} # sample_id_2_max_temp for maximum daily temperature using Sample ID as keys
    sample_id_2_min_temp = {} # sample_id_2_min_temp for minimum daily temperature using Sample ID as keys
    sample_id_2_mean_temp = {}  # sample_id_2_mean_temp for average daily temperature using Sample ID as keys
    sample_id_2_percipitation = {}  # sample_id_2_percipitation for daily total percipitation using Sample ID as keys
    sample_id_2_wind_dir = {}   # sample_id_2_wind_dir for the wind direction using Sample ID as keys
    sample_id_2_max_wind_spd = {}  # sample_id_2_max_wind_spd for daily maximum wind speed using Sample ID as keys
    aTester = {}

    #Begin Reading CSV File
    with open(meta_data, "rb") as meta_data_file:
        meta_data_reader = csv.reader(meta_data_file, delimiter=",")
        meta_data_reader.next() #Skip the Header
        #Iterate over the metaData file, read it row by row, and store the information to the appropriate dictionaries using either sff_mid as keys or Sample ID as keys
        for row in meta_data_reader:
            sff_mid = str(row[3]) + "." + 'MID' + str(row[6])
            sff_mid_2_sample_id[sff_mid] = row[7] #Generate the sff_mid number using row[3] and row[6]
            sample_id_2_collector[row[7]] = row[8]
            sample_id_2_year[row[7]] = row[9]
            sample_id_2_week[row[7]] = row[10]
            sample_id_2_city[row[7]] = row[13]
            sample_id_2_province[row[7]] = row[14]
            aTester.setdefault(row[7], []).append(row[8])
            aTester.setdefault(row[7], []).append(row[9])
            aTester.setdefault(row[7], []).append(row[10])
            aTester.setdefault(row[7], []).append(row[13])
            aTester.setdefault(row[7], []).append(row[14])
            try:
                sample_id_2_dates[row[7]] = str(row[33]) + '-' + str(row[9])
            except IndexError:
                pass
            try:
                sample_id_2_max_temp[row[7]] = row[34]
            except IndexError:
                pass
            try:
                sample_id_2_min_temp[row[7]] = row[35]
            except IndexError:
                pass
            try:
                sample_id_2_mean_temp[row[7]] = row[36]
            except IndexError:
                pass
            try:
                sample_id_2_percipitation[row[7]] = row[41]
            except IndexError:
                pass
            try:
                sample_id_2_wind_dir[row[7]] = row[43]
            except IndexError:
                pass
            try:
                sample_id_2_max_wind_spd[row[7]] = row[44]
            except IndexError:
                pass
            
        
    return sff_mid_2_sample_id, sample_id_2_collector, sample_id_2_year, sample_id_2_week, sample_id_2_city, sample_id_2_province, sample_id_2_dates, sample_id_2_max_temp, sample_id_2_min_temp, sample_id_2_mean_temp, sample_id_2_percipitation, sample_id_2_wind_dir, sample_id_2_max_wind_spd

#The summary_log_reader function reads alll the summary log files row by row,
# When it finds a row with the matching rank and name, it stores the name of the SFF.MID number to a list
# Summary Log Files, String, String --> List
def summary_log_reader(summary_logs, name, rank):
    
    #summary_logs = summary_logs + '*.summarylog'
    list_Of_summary_logs = glob.glob(summary_logs) # Using glob to read multiple files in a directory
    all_summary_logs = [] #This list store 

    # Read each file of the summary log
    for aFile in list_Of_summary_logs:
        file_name = str(aFile)
        print 'Reading', file_name #Tells the user which file it is reading right now.
        sep = '.'
        mid_num = file_name.split(sep)[2]
        mid_num = str(mid_num).strip()   #Extract the SFF.MID number from the file name
        with open(aFile, "rb") as a_summary_log:
            aReader = csv.reader(a_summary_log, delimiter = "\t")
            # Check if each row in a file has a matching name and rank with the input from the user
            for row in aReader:
                if rank == 'kingdom':
                    try:
                        if row[13] == name:
                            aName = str(row[0]).strip() + '.' + mid_num
                            all_summary_logs.append(aName)
                            print 'Found', row[13]
                            print 'Found', len(all_summary_logs), 'sequences'
                        else:
                            pass
                    except IndexError:
                        pass
                if rank == 'phylum':
                    try:
                        if row[14] == name:
                            aName = str(row[0]).strip() + '.' + mid_num
                            all_summary_logs.append(aName)
                            print 'Found', row[14]
                            print 'Found', len(all_summary_logs), 'sequences'
                        else:
                            pass
                    except IndexError:
                        pass
                if rank == 'class':
                    try:
                        if row[15] == name:
                            aName = str(row[0]).strip() + '.' + mid_num
                            all_summary_logs.append(aName)
                            print 'Found', row[15]
                            print 'Found', len(all_summary_logs), 'sequences'
                        else:
                            pass
                    except IndexError:
                        pass
                if rank == 'order':
                    try:
                        if row[16] == name:
                            aName = str(row[0]).strip() + '.' + mid_num
                            all_summary_logs.append(aName)
                            print 'Found', row[16]
                            print 'Found', len(all_summary_logs), 'sequences'
                        else:
                            pass
                    except IndexError:
                        pass
                if rank == 'family':
                    try:
                        if row[17] == name:
                            aName = str(row[0]).strip() + '.' + mid_num
                            all_summary_logs.append(aName)
                            print 'Found', row[17]
                            print 'Found', len(all_summary_logs), 'sequences'
                        else:
                            pass
                    except IndexError:
                        pass
                if rank == 'genus':
                    try:
                        if row[18] == name:
                            aName = str(row[0]).strip() + '.' + mid_num
                            all_summary_logs.append(aName)
                            print 'Found', row[18]
                            print 'Found', len(all_summary_logs), 'sequences'
                        else:
                            pass
                    except IndexError:
                        pass
                if rank == 'species':
                    try:
                        sep = ' '
                        genus = (row[19].split(sep)[0]).strip()
                        species = (row[19].split(sep)[1]).strip()
                        the_name_to_check = genus + ' ' + species
                        if name == the_name_to_check:
                            aName = str(row[0]).strip() + '.' + mid_num
                            all_summary_logs.append(aName)
                            print 'Found', the_name_to_check
                            print 'Found', len(all_summary_logs), 'sequences'
                        else:
                            pass
                    except IndexError:
                        pass
    return all_summary_logs 
        
# The ofp_reader function reads the output from the Oligo Fishing Pipeline (OFP) and stores the Sample ID number of each sequence to a list
# Make use of BioPython's feature
# Text in fasta format --> List
def ofp_reader(ofrFile):
    print 'Reading Oligo Fishing Experimenta Data...'
    seq_record = SeqIO.parse(ofrFile, "fasta") #Using Biopython module to parse the fasta file

    list_of_sample_ids = []

    # Read each fasta sequence
    for i in seq_record:
        descr = i.description
        a_sequence_property = dict(x.split(":") for x in re.findall(r"\w+:\s\w+", descr))
        # Each fasta sequences has multiple dscriptions. This breaks down the description into multiple categories using regular expression.
        # The categories are broken up into dictionaries temporarily 
        if a_sequence_property['Project'] == ' SAGES':
            list_of_sample_ids.append(a_sequence_property['Sample'])
        else:
            continue

    list_of_sample_ids = [i.strip() for i in list_of_sample_ids] #Remove empty spaces from the Sample IDs
    
    return list_of_sample_ids

# This function removes keys of a dictionary that are not present in the lists
# Dict, List, List -> Dict
def rm_redundant_meta(aDict, aList, aList2):

    aSet = set(aList)
    aSet2 = set(aList2)
    
    for k in aDict.keys():
        if (k not in aSet) and (k not in aSet2):
            del aDict[k]

    return aDict
    
# This is the main function that calls the other functions
# Make sure of the pandas module
# It produces a CSV file that compare the oligo fishing and automated identification result.
# It takes the rank, name, the summary logs, the metaData file, the OFP output
# String, string, summary log files, metaData File, OFP fasta text file --> CSV file or DataFrame
def analyzer(rank, name, path_to_summary_logs, path_to_metaData, path_to_ofr_result, rm_redundant):

    startTime = datetime.now()
    print 'Now opening all the summary logs and store to memory'
    #Call summary_log_reader function to read all the summary logs using the rank and name
    summary_logs = summary_log_reader(path_to_summary_logs, name, rank)
    #Call oligoFishingReader to read the results from the oligo fishing experiment data
    list_of_orf_elements = ofp_reader(path_to_ofr_result) 
    multi_dicts = meta_reader(path_to_metaData) 
    meta_sample_ids = multi_dicts[0]
    meta_collectors = multi_dicts[1]
    meta_years = multi_dicts[2]
    meta_weeks = multi_dicts[3]
    metaCities = multi_dicts[4]
    meta_provinces = multi_dicts[5]
    meta_dates = multi_dicts[6]
    meta_max_temp = multi_dicts[7]
    meta_min_temp = multi_dicts[8]
    meta_mean_temp = multi_dicts[9]
    meta_percipitation = multi_dicts[10]
    meta_wind_dir = multi_dicts[11]
    meta_max_wind_spd = multi_dicts[12]

    # Convert the list of sff.mid number from the summary_log_reader function to a list of Sample IDs
    listOfAutoElements = []
    for i in summary_logs:
        try: 
            aSampleID = meta_sample_ids[i]
            listOfAutoElements.append(aSampleID)   
        except KeyError:
            pass

    # Control if the Samples that don't have the indicates organism should be included in the final output data sheet
    if rm_redundant == 'Yes':
        meta_collectors = rm_redundant_meta(meta_collectors, listOfAutoElements, list_of_orf_elements)
        meta_years = rm_redundant_meta(meta_years, listOfAutoElements, list_of_orf_elements)
        meta_weeks = rm_redundant_meta(meta_weeks, listOfAutoElements, list_of_orf_elements)
        metaCities = rm_redundant_meta(metaCities, listOfAutoElements, list_of_orf_elements)
        meta_provinces = rm_redundant_meta(meta_provinces, listOfAutoElements, list_of_orf_elements)
        meta_dates = rm_redundant_meta(meta_dates, listOfAutoElements, list_of_orf_elements)
        meta_max_temp = rm_redundant_meta(meta_max_temp, listOfAutoElements, list_of_orf_elements)
        meta_min_temp = rm_redundant_meta(meta_min_temp, listOfAutoElements, list_of_orf_elements)
        meta_mean_temp = rm_redundant_meta(meta_mean_temp, listOfAutoElements, list_of_orf_elements)
        meta_percipitation = rm_redundant_meta(meta_percipitation, listOfAutoElements, list_of_orf_elements)
        meta_wind_dir = rm_redundant_meta(meta_wind_dir, listOfAutoElements, list_of_orf_elements)
        meta_max_wind_spd = rm_redundant_meta(meta_max_wind_spd, listOfAutoElements, list_of_orf_elements)
    else:
        pass

    # Count the appearance of each element in the list of samples IDs and output the results into dictionaries
    num_of_ofp_elements = Counter(list_of_orf_elements)
    num_of_auto_elements = Counter(listOfAutoElements)
    
    # Converting all the data into a dataframe using pandas
    print 'Generating DataFrame and CSV output...'
    temp = {'Provinces': pd.Series(meta_provinces), 'Cities': pd.Series(metaCities), 'Collectors': pd.Series(meta_collectors), 'Years': pd.Series(meta_years), 'Weeks': pd.Series(meta_weeks), 'Dates': pd.Series(meta_dates), 'Max Temperature': pd.Series(meta_max_temp), 'Min Temperature': pd.Series(meta_min_temp), 'Mean Temperature': pd.Series(meta_mean_temp), 'Total Percipitation': pd.Series(meta_percipitation), 'Wind Direction': pd.Series(meta_wind_dir), 'Max Wind Speed': pd.Series(meta_max_wind_spd), 'Auto ID': pd.Series(num_of_auto_elements), 'OFP Result': pd.Series(num_of_ofp_elements)}
    df = pd.DataFrame(temp, columns=['Provinces', 'Cities', 'Collectors', 'Years', 'Weeks', 'Dates', 'Max Temperature', 'Min Temperature', 'Mean Temperature', 'Total Percipitation', 'Wind Direction', 'Max Wind Speed', 'Auto ID', 'OFP Result'])
    df = df.fillna(value=0)
    df['Log(Auto ID + 1)'] = df['Auto ID'].apply(np.log1p)
    df['Log(OFP Result + 1)'] = df['OFP Result'].apply(np.log1p)
    file_name = rank + '_' + name + '_' + 'Auto vs OFP.csv'
    df.to_csv(file_name) #Output: Writing the dataframe to a CSV file
    print 'Took', (datetime.now()-startTime)
    return df
        
    

if __name__=='__main__':
##    analyzer(sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    item1 = analyzer('species', 'Staphylococcus aureus', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Staphylococcus_aureus.Staphylococcus_aureus.txt.withmeta.fna.fasta', 'No')
##    item2 = analyzer('species', 'Bacillus anthracis', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Bacillus_anthracis_cand1.Bacillus_anthracis_cand1.txt.withmeta.fna.fasta', 'No')
##    item3 = analyzer('species', 'Bacillus herbersteinensis', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Bacillus_anthracis_cand1.Bacillus_anthracis_cand1.txt.withmeta.fna.fasta', 'No')
##    item4 = analyzer('species', 'Staphylococcus epidermidis', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Staphylococcus_epidermidis.Staphylococcus_epidermidis.txt.withmeta.fna.fasta', 'No')
##    item5 = analyzer('species', 'Clostridium botulinum', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Clostridium_botulinum.Clostridium_botulinum.txt.withmeta.fna.fasta', 'No')
##    item6 = analyzer('species', 'Clostridium perfringens', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Clostridium_perfringens.Clostridium_perfringens.txt.withmeta.fna.fasta', 'No')
##    item7 = analyzer('phylum', 'Firmicutes', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Firmicutes_oligo1.Firmicutes.txt.withmeta.fna.fasta', 'No')
##    item8 = analyzer('genus', 'Mycobacterium', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Mycobacterium_spp.Mycobacterium_spp.txt.withmeta.fna.fasta', 'No')
##    item9 = analyzer('genus', 'Streptococcus', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Streptococcus_spp.Streptococcus_spp.txt.withmeta.fna.fasta', 'No')
##    item10 = analyzer('genus', 'Bacillus', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Bacillus_spp.Bacillus_spp.txt.withmeta.fna.fasta', 'No')
##    item11 = analyzer('genus', 'Staphylococcus', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Staphylococcus_spp.Staphylococcus_spp.txt.withmeta.fna.fasta', 'No')
##    item12 = analyzer('genus', 'Rickettsia', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Rickettsia_spp.Rickettsia_spp.txt.withmeta.fna.fasta', 'No')
##    item13 = analyzer('genus', 'Clostridium', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Clostridium_spp.Clostridium_spp.txt.withmeta.fna.fasta', 'No')




