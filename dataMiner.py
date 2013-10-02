# By Stephen Fung
# stephen.fung@agr.gc.ca
# The script reads the Oligo Fishing Pipeline (OFP) result and Automated Identification (Auto ID) result based on a given rank and a name.
# It outputs the Auto ID and OFP data and return the climate, geography...etc in a csv file that can be opened in MS Excel by the biologists
# The script requires numpy, biopython and pandas to be installed with python 2.7
#Testing 2013 Fall

import csv
import re
import sys
import glob
import os
import math
from collections import Counter
from collections import defaultdict
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO



# The metaReader function reads the metaData File and return multiple dictionaries:
# The dictionaries uses Sample ID as keys for the relevant information
# CSV File --> Dictionaries
def meta_reader(meta_data):

    sff_mid_2_sample_id = {} # sff_mid_2_Sample for finding out Sample ID using sff_mid as keys
    sample_ids_2_meta_info = defaultdict(list) # a hashtable using the ID of the samples as key to store: collector, year, week, locations, and climate information as value (in a list)

    #Begin Reading CSV File
    with open(meta_data, "rb") as meta_data_file:
        meta_data_reader = csv.reader(meta_data_file, delimiter=",")
        meta_data_reader.next() #Skip the Header
        #Iterate over the metaData file, read it row by row, and store the information to the appropriate dictionaries using either sff_mid as keys or Sample ID as keys
        for row in meta_data_reader:
            sff_mid = str(row[3]) + "." + 'MID' + str(row[6])
            sff_mid_2_sample_id[sff_mid] = row[7] #Generate the sff_mid number using row[3] and row[6]

            a_sample_id = row[7]  #Sample ID of each sample
            a_collector = row[8]    #The collector used for the sample
            a_year = row[9] # The year the sample was collected
            a_week = row[10] # The week number the sample was collected
            a_city = row[13]    # The city the sample was collected
            a_province = row[14] # The province the sample was collected
            try:
                a_date = str(row[33]) + '-' + str(row[9])  #Using row[33] and row[9] to construct the actual date of the sample
            except IndexError:  # Some cells can be empty
                a_date = None
            try:
                a_max_temp = row[34]    # Maximum Temperature
            except IndexError:
                a_max_temp = None
            try:
                a_min_temp = row[35]    # Minimum Temperature
            except IndexError:
                a_min_temp = None
            try:
                a_mean_temp = row[36]   # Mean Temperature
            except IndexError:
                a_mean_temp = None
            try:
                a_percipitation = row[41]   # Percipitation
            except IndexError:
                a_percipitation = None
            try:
                a_wind_dir = row[43]    # Wind Direction
            except IndexError:
                a_wind_dir = None
            try:
                a_max_wind_spd = row[44]    # Maximum Wind Speed of the Date
            except IndexError:
                a_max_wind_spd = None

            # Group the information to a list and add it to the defaultdict
            a_row = [a_collector, a_year, a_week, a_city, a_province, a_date, a_max_temp, a_min_temp, a_mean_temp, a_percipitation, a_wind_dir, a_max_wind_spd]
            
            if a_sample_id not in sample_ids_2_meta_info:
                sample_ids_2_meta_info[a_sample_id].extend(a_row) 
                        
        
    return sff_mid_2_sample_id, sample_ids_2_meta_info

#The summary_log_reader function reads alll the summary log files row by row,
# When it finds a row with the matching rank and name, it stores the name of the SFF.MID number to a list
# Summary Log Files, String, String --> List
def summary_log_reader(summary_logs, name, rank):
    
    #summary_logs = summary_logs + '*.summarylog'
    list_Of_summary_logs = glob.glob(summary_logs) # Using glob to read multiple files in a directory
    all_summary_logs = [] #This list store 

    # Read each file of the summary log
    for a_file in list_Of_summary_logs:
        file_name = str(a_file)
        print 'Reading', file_name #Tells the user which file it is reading right now.
        sep = '.'
        mid_num = file_name.split(sep)[2]
        mid_num = str(mid_num).strip()   #Extract the mid number from the file name
        with open(a_file, "rb") as a_summary_log:
            aReader = csv.reader(a_summary_log, delimiter = "\t")
            # Check if each row in a file has a matching name and rank with the input from the user
            # Store the sff mid number to all_summary_logs each time the name of the organism appears in the row
            for row in aReader:
                if rank == 'kingdom':
                    try:
                        if row[13] == name:
                            aName = str(row[0]).strip() + '.' + mid_num # Generate the sff mid number here
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
def analyzer(rank, name, path_to_summary_logs, path_to_meta_data, path_to_ofr_result, rm_redundant):

    startTime = datetime.now()
    print 'Now opening all the summary logs and store to memory'
    #Call summary_log_reader function to read all the summary logs using the rank and name
    summary_logs = summary_log_reader(path_to_summary_logs, name, rank)
    #Call oligoFishingReader to read the results from the oligo fishing experiment data
    list_of_orf_elements = ofp_reader(path_to_ofr_result) 
    meta_sample_ids, sample_ids_2_meta_info = meta_reader(path_to_meta_data) 


    # Convert the list of sff mid number from the summary_log_reader function to a list of Sample IDs
    list_of_auto_elements = []
    for i in summary_logs:
        try: 
            a_sample_id = meta_sample_ids[i]
            list_of_auto_elements.append(a_sample_id)   
        except KeyError:
            pass

    # Control if the Samples that don't have the indicates organism should be included in the final output data sheet
    if rm_redundant == 'Yes':
        sample_ids_2_meta_info = rm_redundant_meta(sample_ids_2_meta_info, list_of_auto_elements, list_of_orf_elements)
    else:
        pass

    # Count the appearance of each element in the list of samples IDs and output the results into dictionaries
    num_of_ofp_elements = Counter(list_of_orf_elements)
    num_of_auto_elements = Counter(list_of_auto_elements)
    
    # Converting all the data into a dataframe using pandas
    print 'Generating DataFrame and CSV output...'
    meta_df = pd.DataFrame(sample_ids_2_meta_info, index=['Collector', 'Year', 'Week', 'Cities', 'Provinces', 'Dates', 'Max Temperature', 'Min Temperature', 'Mean Temperature', 'Total Percipitation', 'Wind Direction', 'Max Wind Speed'])
    meta_df = meta_df.transpose()
    result_dict = {'Auto ID': pd.Series(num_of_auto_elements), 'OFP Result': pd.Series(num_of_ofp_elements)}
    result_df = pd.DataFrame(result_dict)
    final_df = meta_df.join(result_df)
    final_df = final_df.fillna(value=0) #Combine the meta dataframe with the result dataframe
    final_df['Log(Auto ID + 1)'] = final_df['Auto ID'].apply(np.log1p) #Apply log10 to the counts
    final_df['Log(OFP Result + 1)'] = final_df['OFP Result'].apply(np.log1p)
    file_name = rank + '_' + name + '_' + 'Auto vs OFP.csv'
    final_df.to_csv(file_name) #Output: Writing the dataframe to a CSV file
    print 'Took', (datetime.now()-startTime)
    return final_df
        
    

if __name__=='__main__':
    analyzer(sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
##    item1 = analyzer('species', 'Staphylococcus aureus', './rs75_Summary_Log/*.summarylog', './input/454_sample_summary_with_climate.csv', './input/Staphylococcus_aureus.Staphylococcus_aureus.txt.withmeta.fna.fasta', 'No')
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




