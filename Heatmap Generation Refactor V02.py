
import csv
import math
import numpy as np 
import pandas
import glob
import pylab as pl
import matplotlib.pyplot as plt


"""
#Test Reading the Data directly...

d = {}


for i in glob.glob('*.tsv'):
    #Process the file title to make it cleaner
    title = i
    title = title[22:]
    title = title.split('.tsv')[0]
    title = title.split('.sff')[0]

    #open each tsv file using the csv reader with tab as the delimiter
    temp = csv.reader(open(i, 'r'), delimiter = "\t")

    #Skipping the first two title lines
    temp.next()
    temp.next()


    for row in temp:
        try:
            indexKey= row[2] + '.' + row[0] + '.' + row[3].strip('txid:')

            if d.has_key(indexKey):
                d[indexKey][title] = row[1]
            else:
                d[indexKey] = {}
                d[indexKey][title] = row[1]

        except IndexError:
            continue
"""


#Pathogens Genus and Species Index
pathogensIndexByGenus = []
pathogensIndexBySpecies = []
rawData = open('Taxonomy Result.csv', 'r')
rawDataReader = csv.reader(rawData, delimiter = ',')
for row in rawDataReader:
    if row[0] == "genus":
        if (row[1] == "Lactobacillus"
            or row[1] == "Escherichia"
            or row[1] == "Clostridium"
            or row[1] == "Campylobacter"
            or row[1] == "Salmonella"
            or row[1] == "Bifidobacterium"
            or row[1] == "Staphylococcus"
            or row[1] == "Bacillus"
            or row[1] == "Bordetella"
            or row[1] == "Brucella"
            or row[1] == "Chlamydia"
            or row[1] == "Chlamydophila"
            or row[1] == "Corynebacterium"
            or row[1] == "Enterococcus"
            or row[1] == "Francisella"
            or row[1] == "Haemophilus"
            or row[1] == "Mycoplasma"
            or row[1] == "Helicobacter"
            or row[1] == "Legionella"
            or row[1] == "Leptospira"
            or row[1] == "Listeria"
            or row[1] == "Mycobacterium"
            or row[1] == "Neisseria"
            or row[1] == "Pseudomonas"
            or row[1] == "Rickettsia"
            or row[1] == "Shigella"
            or row[1] == "Streptococcus"
            or row[1] == "Treponema"
            or row[1] == "Vibrio"
            or row[1] == "Yersinia"
            or row[1] == "Brucella"):
            pathogensIndexByGenus.append(row[1])
            
    if (row[1] == 'Bacillus anthracis'
        or row[1] == 'Bacillus cereus'
        or row[1] == 'Bacillus thuringiensis'
        or row[1] == 'Bordetella pertussis'
        or row[1] == 'Brucella abortus'
        or row[1] == 'Brucella canis'
        or row[1] == 'Brucella melitensis'
        or row[1] == 'Brucella suis'
        or row[1] == 'Campylobacter jejuni'
        or row[1] == 'Coxiella burnetii'
        or row[1] == 'Francisella tularensis'
        or row[1] == 'Mycobacterium tuberculosis complex'
        or row[1] == 'Clostridium botulinum'
        or row[1] == 'Clostridium perfringens'
        or row[1] == '[Clostridium] difficile'
        or row[1] == 'Clostridium tetani'
        or row[1] == 'Staphylococcus aureus'
        or row[1] == 'Bacillus subtilis'
        or row[1] == 'Bacillus subtilis group'
        or row[1] == 'Bacillus prodigiosus'
        or row[1] == 'Bacillus patchiness'
        or row[1] == 'Bacillus violaceus'
        or row[1] == 'Chlamydia pneumoniae'
        or row[1] == 'Chlamydia trachomatis'
        or row[1] == 'Chlamydophila psittaci'
        or row[1] == 'Cytophaga allerginae'
        or row[1] == 'Corynebacterium diphtheriae'
        or row[1] == 'Enterobacter cloacae'
        or row[1] == 'Enterobacter cloacae complex'
        or row[1] == 'Erwinia herbicola'
        or row[1] == 'Escherichia coli'
        or row[1] == 'Klebisella planticola'
        or row[1] == 'Legionella pneumophila'
        or row[1] == 'Mycoplasma pneumonaie'
        or row[1] == 'Mycobacteria avium'
        or row[1] == 'Mycobacterium leprae'
        or row[1] == 'Mycobacterium ulcerans'
        or row[1] == 'Neisseria gonorrhoeae'
        or row[1] == 'Neisseria_meningitidis'
        or row[1] == 'Pseudomonas aeruginosa group'
        or row[1] == 'Rickettsia rickettsii'
        or row[1] == 'Salmonella typhi'
        or row[1] == 'Salmonella typhimurium'
        or row[1] == 'Shigella sonnei'
        or row[1] == 'Treponema pallidum'
        or row[1] == 'Francisella tularensis'
        or row[1] == 'Serratia marcescens'
        or row[1] == 'Staphyococcus albus'
        or row[1] == 'Staphylococcus epidermidis'
        or row[1] == 'Staphylococcus saprophyticus'
        or row[1] == 'Streptococcus agalactiae'
        or row[1] == 'Streptococcus pneumoniae'
        or row[1] == 'Streptococcus pyogenes'
        or row[1] == "Yersinia pestis"
        or row[1] == 'Vibrio cholerae'
        or row[1] == 'Shigella dysenteriae'
        or row[1] == 'Burkholderia mallei'
        or row[1] == 'Haemophilus influenzae'
        or row[1] == 'Helicobacter pylori'
        or row[1] == 'Leptospira interrogans'
        or row[1] == 'Listeria monocytogenes'):
            pathogensIndexBySpecies.append(row[1])
    
        
   

#Reading all the taxonomy into a data frame
            
taxonomyDataFrame = pandas.io.parsers.read_csv("Taxonomy Result.csv", index_col=[0,1,2], header=5)


# dataframe, string, listofstring, listofstring -> dataframe
def dataFrameSlicer(df, x, lopg=pathogensIndexByGenus, lops=pathogensIndexBySpecies):
    if x == 'species':
        dfSpecies = df.ix['species']
        dfSpeciesPathogens = dfSpecies.ix[lops]    #Only the Species in the lops is used
        return dfSpeciesPathogens
    if x == 'genus':
        dfGenus =  df.ix['genus']
        dfPathogenGenus = dfGenus.ix[lopg]   #Only the Genus in the lopg is used
        dfPathogenGenus = dfPathogenGenus.drop(('Bacillus',55087)) #Suppose to be some kind of Eukaryotes...showing up in bacteria for some reason...
        return dfPathogenGenus
    if x == 'phylum':
        dfPhylum = df.ix['phylum']
        dfSuperPhylum = df.ix['superphylum']
        dfSubPhylum = df.ix['subphylum']
        dfPhylum = dfPhylum.append(dfSuperPhylum)
        dfPhylum = dfPhylum.append(dfSubPhylum)
        return dfPhylum
    else:
        print 'You have to enter a proper data frame and parameter to the function!'


#Phylum DataFrame
taxonomyPhylumDataFrame = dataFrameSlicer(taxonomyDataFrame, 'phylum')

#Pathogens (Genus) DataFrame
PathogensGenusDataFrame = dataFrameSlicer(taxonomyDataFrame, 'genus')

#Pathogens (Species) DataFrame
PathogensSpeciesDataFrame = dataFrameSlicer(taxonomyDataFrame, 'species')




#Data Frame of Sample Information
indexDataFrame = pandas.io.parsers.read_csv("Index V23.csv")
indexDataFrame = indexDataFrame[(indexDataFrame.Project == 'SAGES')]
indexDataFrame = indexDataFrame[(indexDataFrame.Target == 'Bacteria')]

#Building a dictionary of Sample ID: location, date, collector, target...
locationsIndex = {}
dateIndex ={}
collectorIndex = {}
targetIndex = {}

for r in indexDataFrame.itertuples():
    sampleID = r[2] + '.MID' + str(int(r[4]))

    if sampleID != 'HI4EL9301.MID49' and sampleID != 'G42MIWK02.MID43':    #Exception: This file is empty in the database for some reason....
        locationsIndex[sampleID] = r[10]
        dateIndex[sampleID] = r[9]
        collectorIndex[sampleID] = r[6]
        targetIndex[sampleID] = r[16]

#Get all the Dates
listOfAllDates = []
for key, value in dateIndex.iteritems():
    if value not in listOfAllDates:
        listOfAllDates.append(value)

listOfAllDates = sorted(listOfAllDates)
listOfAllDates = listOfAllDates[1:]



# dataframe, string, string, listofstring, listofstring, listofstring, listofstring -> dataframe

def transformDataFrameWithRequest(df, location, collector, dates=listOfAllDates, id1=locationsIndex, id2=dateIndex, id3=collectorIndex):
    collectorSet = set()
    for key, value in id3.iteritems():
        if value == collector:
            collectorSet.add(key)
    
    # Build location first according to location inputed by the user
    locationList = []
    for key, value in id1.iteritems():
        if value == location:
            locationList.append(key)
    # Cross check the locations with collector requsted by the user
    collectorAndLocation = [i for i in locationList if i in collectorSet]
    
    timeList = []
    for i in collectorAndLocation:
        if i in id2.keys():
            timeList.append(id2[i])
    
    dfUpdatedWithLocationAndCollector = df[collectorAndLocation]
    dfUpdatedWithLocationAndCollector.columns = timeList
    dfUpdatedWithLocationAndCollector = dfUpdatedWithLocationAndCollector.fillna(value=0)
    dfUpdatedWithLocationAndCollector = dfUpdatedWithLocationAndCollector.groupby(axis=1, level=0).sum()
    #dfUpdatedWithLocationAndCollector = dfUpdatedWithLocationAndCollector.reindex(columns=dates)
    dfUpdatedWithLocationAndCollector = dfUpdatedWithLocationAndCollector.sort_index(axis=1)
    dfUpdatedWithLocationAndCollector = dfUpdatedWithLocationAndCollector.sort_index(axis=0)
    dfUpdatedWithLocationAndCollector = dfUpdatedWithLocationAndCollector.apply(np.log1p)        
    
    return dfUpdatedWithLocationAndCollector
        

# Phylum Data for Graphs
taxonomyPhylumDataFrame_AB_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Lacombe_AB', 'JB')
taxonomyPhylumDataFrame_AB_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Lacombe_AB', 'YE')
taxonomyPhylumDataFrame_AB_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Lacombe_AB', 'BK')

taxonomyPhylumDataFrame_BC_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Agassiz_BC', 'JB')
taxonomyPhylumDataFrame_BC_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Agassiz_BC', 'YE')
taxonomyPhylumDataFrame_BC_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Agassiz_BC', 'BK')

taxonomyPhylumDataFrame_SK_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Saskatoon_SK', 'JB')
taxonomyPhylumDataFrame_SK_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Saskatoon_SK', 'YE')
taxonomyPhylumDataFrame_SK_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Saskatoon_SK', 'BK')

taxonomyPhylumDataFrame_MB_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Morden_MB', 'JB')
taxonomyPhylumDataFrame_MB_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Morden_MB', 'YE')
taxonomyPhylumDataFrame_MB_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Morden_MB', 'BK')

taxonomyPhylumDataFrame_Ottawa_ON_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Ottawa_ON', 'JB')
taxonomyPhylumDataFrame_Ottawa_ON_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Ottawa_ON', 'YE')
taxonomyPhylumDataFrame_Ottawa_ON_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Ottawa_ON', 'BK')

taxonomyPhylumDataFrame_Ridgetown_ON_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Ridgetown_ON', 'JB')
taxonomyPhylumDataFrame_Ridgetown_ON_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Ridgetown_ON', 'YE')
taxonomyPhylumDataFrame_Ridgetown_ON_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Ridgetown_ON', 'BK')

taxonomyPhylumDataFrame_Woodslee_ON_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Woodslee_ON', 'JB')
taxonomyPhylumDataFrame_Woodslee_ON_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Woodslee_ON', 'YE')
taxonomyPhylumDataFrame_Woodslee_ON_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Woodslee_ON', 'BK')

taxonomyPhylumDataFrame_QC_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'SainteClotilde_QC', 'JB')
taxonomyPhylumDataFrame_QC_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'SainteClotilde_QC', 'YE')
taxonomyPhylumDataFrame_QC_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'SainteClotilde_QC', 'BK')

taxonomyPhylumDataFrame_PEI_JB = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Harrington_PEI', 'JB')
taxonomyPhylumDataFrame_PEI_YE = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Harrington_PEI', 'YE')
taxonomyPhylumDataFrame_PEI_BK = transformDataFrameWithRequest(taxonomyPhylumDataFrame, 'Harrington_PEI', 'BK')




#Pathogens Genus Data for Graph
AB_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Lacombe_AB', 'JB')
AB_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Lacombe_AB', 'YE')
AB_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Lacombe_AB', 'BK')

BC_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Agassiz_BC', 'JB')
BC_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Agassiz_BC', 'YE')
BC_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Agassiz_BC', 'BK')

SK_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Saskatoon_SK', 'JB')
SK_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Saskatoon_SK', 'YE')
SK_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Saskatoon_SK', 'BK')

MB_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Morden_MB', 'JB')
MB_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Morden_MB', 'YE')
MB_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Morden_MB', 'BK')

Ottawa_ON_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Ottawa_ON', 'JB')
Ottawa_ON_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Ottawa_ON', 'YE')
Ottawa_ON_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Ottawa_ON', 'BK')

Ridgetown_ON_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Ridgetown_ON', 'JB')
Ridgetown_ON_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Ridgetown_ON', 'YE')
Ridgetown_ON_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Ridgetown_ON', 'BK')

Woodslee_ON_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Woodslee_ON', 'JB')
Woodslee_ON_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Woodslee_ON', 'YE')
Woodslee_ON_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Woodslee_ON', 'BK')


QC_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'SainteClotilde_QC', 'JB')
QC_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'SainteClotilde_QC', 'YE')
QC_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'SainteClotilde_QC', 'BK')


PEI_JB_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Harrington_PEI', 'JB')
PEI_YE_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Harrington_PEI', 'YE')
PEI_BK_Pathogens_Genus_DataFrame = transformDataFrameWithRequest(PathogensGenusDataFrame, 'Harrington_PEI', 'BK')



#Pathogens Species Data for Graph

AB_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Lacombe_AB', 'JB')
AB_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Lacombe_AB', 'YE')
AB_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Lacombe_AB', 'BK')


BC_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Agassiz_BC', 'JB')
BC_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Agassiz_BC', 'YE')
BC_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Agassiz_BC', 'BK')


SK_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Saskatoon_SK', 'JB')
SK_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Saskatoon_SK', 'YE')
SK_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Saskatoon_SK', 'BK')


MB_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Morden_MB', 'JB')
MB_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Morden_MB', 'YE')
MB_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Morden_MB', 'BK')


Ottawa_ON_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Ottawa_ON', 'JB')
Ottawa_ON_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Ottawa_ON', 'YE')
Ottawa_ON_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Ottawa_ON', 'BK')


Ridgetown_ON_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Ridgetown_ON', 'JB')
Ridgetown_ON_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Ridgetown_ON', 'YE')
Ridgetown_ON_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Ridgetown_ON', 'BK')


Woodslee_ON_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Woodslee_ON', 'JB')
Woodslee_ON_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Woodslee_ON', 'YE')
Woodslee_ON_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Woodslee_ON', 'BK')


QC_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'SainteClotilde_QC', 'JB')
QC_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'SainteClotilde_QC', 'YE')
QC_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'SainteClotilde_QC', 'BK')


PEI_JB_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Harrington_PEI', 'JB')
PEI_YE_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Harrington_PEI', 'YE')
PEI_BK_Pathogens_Species_DataFrame = transformDataFrameWithRequest(PathogensSpeciesDataFrame, 'Harrington_PEI', 'BK')




# dataframe, string, string, string, string, string -> graph

def drawHeatMap(df, city, province, collector, classtype, color, titleposy):
    try:
        thePlot = pl.matshow(df.values, cmap='PuBuGn')
        pl.colorbar(thePlot, orientation='vertical')
        aTitle = classtype + ' Composition Changes Over Time in ' + city + ', ' + province + '\n' + collector + ' collector. ' + 'rs100'
        pl.title(aTitle, x=0.5, y=titleposy, style='oblique', weight='bold')
        pl.xlabel('Collection Time')
        #pl.ylabel(classtype + ' and Taxonomy ID')
        pl.xticks(range(len(df.columns)), df.columns, rotation=90)
        pl.yticks(range(len(df.index)), df.index)
    except ZeroDivisionError:
        errorMessage = 'No Data Avaiable for ' + city + ', ' + province + ' with ' + collector + ' collector.'
        print errorMessage
    
        
"""
print 'Processing JB Collector Data...'

drawHeatMap(taxonomyPhylumDataFrame_AB_JB, 'Lacombe', 'AB', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_BC_JB, 'Agassiz', 'BC', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_SK_JB, 'Saskatoon', 'SK', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_MB_JB, 'Morden', 'MB', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Ottawa_ON_JB, 'Ottawa', 'ON', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Ridgetown_ON_JB, 'Ridgetown', 'ON', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Woodslee_ON_JB, 'Woodslee', 'ON', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_QC_JB, 'Sainte Clotilde', 'QC', 'JB', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_PEI_JB, 'Harrington', 'PEI', 'JB', 'Phylum', 'YlGn', 1.1)


drawHeatMap(AB_JB_Pathogens_Genus_DataFrame, 'Lacombe', 'AB', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(BC_JB_Pathogens_Genus_DataFrame, 'Agassiz', 'BC', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(SK_JB_Pathogens_Genus_DataFrame, 'Saskatoon', 'SK', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(MB_JB_Pathogens_Genus_DataFrame, 'Morden', 'MB', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(Ottawa_ON_JB_Pathogens_Genus_DataFrame, 'Ottawa', 'ON', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(Ridgetown_ON_JB_Pathogens_Genus_DataFrame, 'Ridgetown', 'ON', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(Woodslee_ON_JB_Pathogens_Genus_DataFrame, 'Woodslee', 'ON', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(QC_JB_Pathogens_Genus_DataFrame, 'Sainte Clotilde', 'QC', 'JB', 'Genus', 'YlGn', 1.1)
drawHeatMap(PEI_JB_Pathogens_Genus_DataFrame, 'Harrington', 'PEI', 'JB', 'Genus', 'YlGn', 1.1)


drawHeatMap(AB_JB_Pathogens_Species_DataFrame, 'Lacombe', 'AB', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(BC_JB_Pathogens_Species_DataFrame, 'Agassiz', 'BC', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(SK_JB_Pathogens_Species_DataFrame, 'Saskatoon', 'SK', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(MB_JB_Pathogens_Species_DataFrame, 'Morden', 'MB', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(Ottawa_ON_JB_Pathogens_Species_DataFrame, 'Ottawa', 'ON', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(Ridgetown_ON_JB_Pathogens_Species_DataFrame, 'Ridgetown', 'ON', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(Woodslee_ON_JB_Pathogens_Species_DataFrame, 'Woodslee', 'ON', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(QC_JB_Pathogens_Species_DataFrame, 'Sainte Clotilde', 'QC', 'JB', 'Species', 'YlGn', 1.1)
drawHeatMap(PEI_JB_Pathogens_Species_DataFrame, 'Harrington', 'PEI', 'JB', 'Species', 'YlGn', 1.1)



print 'Processing YE Collector Data...'

drawHeatMap(taxonomyPhylumDataFrame_AB_YE, 'Lacombe', 'AB', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_BC_YE, 'Agassiz', 'BC', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_SK_YE, 'Saskatoon', 'SK', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_MB_YE, 'Morden', 'MB', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Ottawa_ON_YE, 'Ottawa', 'ON', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Ridgetown_ON_YE, 'Ridgetown', 'ON', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Woodslee_ON_YE, 'Woodslee', 'ON', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_QC_YE, 'Sainte Clotilde', 'QC', 'YE', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_PEI_YE, 'Harrington', 'PEI', 'YE', 'Phylum', 'YlGn', 1.1)



drawHeatMap(AB_YE_Pathogens_Genus_DataFrame, 'Lacombe', 'AB', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(BC_YE_Pathogens_Genus_DataFrame, 'Agassiz', 'BC', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(SK_YE_Pathogens_Genus_DataFrame, 'Saskatoon', 'SK', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(MB_YE_Pathogens_Genus_DataFrame, 'Morden', 'MB', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(Ottawa_ON_YE_Pathogens_Genus_DataFrame, 'Ottawa', 'ON', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(Ridgetown_ON_YE_Pathogens_Genus_DataFrame, 'Ridgetown', 'ON', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(Woodslee_ON_YE_Pathogens_Genus_DataFrame, 'Woodslee', 'ON', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(QC_YE_Pathogens_Genus_DataFrame, 'Sainte Clotilde', 'QC', 'YE', 'Genus', 'YlGn', 1.1)
drawHeatMap(PEI_YE_Pathogens_Genus_DataFrame, 'Harrington', 'PEI', 'YE', 'Genus', 'YlGn', 1.1)



drawHeatMap(AB_YE_Pathogens_Species_DataFrame, 'Lacombe', 'AB', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(BC_YE_Pathogens_Species_DataFrame, 'Agassiz', 'BC', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(SK_YE_Pathogens_Species_DataFrame, 'Saskatoon', 'SK', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(MB_YE_Pathogens_Species_DataFrame, 'Morden', 'MB', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(Ottawa_ON_YE_Pathogens_Species_DataFrame, 'Ottawa', 'ON', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(Ridgetown_ON_YE_Pathogens_Species_DataFrame, 'Ridgetown', 'ON', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(Woodslee_ON_YE_Pathogens_Species_DataFrame, 'Woodslee', 'ON', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(QC_YE_Pathogens_Species_DataFrame, 'Sainte Clotilde', 'QC', 'YE', 'Species', 'YlGn', 1.1)
drawHeatMap(PEI_YE_Pathogens_Species_DataFrame, 'Harrington', 'PEI', 'YE', 'Species', 'YlGn', 1.1)



print 'Processing BK Collector Data...'


drawHeatMap(taxonomyPhylumDataFrame_AB_BK, 'Lacombe', 'AB', 'BK', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_BC_BK, 'Agassiz', 'BC', 'BK', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_SK_BK, 'Saskatoon', 'SK', 'BK', 'Phylum', 'YlGn',1.1)
drawHeatMap(taxonomyPhylumDataFrame_MB_BK, 'Morden', 'MB', 'BK', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Ottawa_ON_BK, 'Ottawa', 'ON', 'BK', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Ridgetown_ON_BK, 'Ridgetown', 'ON', 'BK', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_Woodslee_ON_BK, 'Woodslee', 'ON', 'BK', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_QC_BK, 'Sainte Clotilde', 'QC', 'BK', 'Phylum', 'YlGn', 1.1)
drawHeatMap(taxonomyPhylumDataFrame_PEI_BK, 'Harrington', 'PEI', 'BK', 'Phylum', 'YlGn', 1.1)



drawHeatMap(AB_BK_Pathogens_Genus_DataFrame, 'Lacombe', 'AB', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(BC_BK_Pathogens_Genus_DataFrame, 'Agassiz', 'BC', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(SK_BK_Pathogens_Genus_DataFrame, 'Saskatoon', 'SK', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(MB_BK_Pathogens_Genus_DataFrame, 'Morden', 'MB', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(Ottawa_ON_BK_Pathogens_Genus_DataFrame, 'Ottawa', 'ON', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(Ridgetown_ON_BK_Pathogens_Genus_DataFrame, 'Ridgetown', 'ON', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(Woodslee_ON_BK_Pathogens_Genus_DataFrame, 'Woodslee', 'ON', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(QC_BK_Pathogens_Genus_DataFrame, 'Sainte Clotilde', 'QC', 'BK', 'Genus', 'YlGn', 1.1)
drawHeatMap(PEI_BK_Pathogens_Genus_DataFrame, 'Harrington', 'PEI', 'BK', 'Genus', 'YlGn', 1.1)


drawHeatMap(AB_BK_Pathogens_Species_DataFrame, 'Lacombe', 'AB', 'BK', 'Species', 'YlGn', 1.1)
drawHeatMap(BC_BK_Pathogens_Species_DataFrame, 'Agassiz', 'BC', 'BK', 'Species', 'YlGn', 1.1)
drawHeatMap(SK_BK_Pathogens_Species_DataFrame, 'Saskatoon', 'SK', 'BK', 'Species', 'YlGn', 1.1)
drawHeatMap(MB_BK_Pathogens_Species_DataFrame, 'Morden', 'MB', 'BK', 'Species', 'YlGn', 1.1)
drawHeatMap(Ottawa_ON_BK_Pathogens_Species_DataFrame, 'Ottawa', 'ON', 'BK', 'Species', 'YlGn', 1.2)
drawHeatMap(Ridgetown_ON_BK_Pathogens_Species_DataFrame, 'Ridgetown', 'ON', 'BK', 'Species', 'YlGn', 1.2)
drawHeatMap(Woodslee_ON_BK_Pathogens_Species_DataFrame, 'Woodslee', 'ON', 'BK', 'Species', 'YlGn', 1.2)
drawHeatMap(QC_BK_Pathogens_Species_DataFrame, 'Sainte Clotilde', 'QC', 'BK', 'Species', 'YlGn', 1.1)
drawHeatMap(PEI_BK_Pathogens_Species_DataFrame, 'Harrington', 'PEI', 'BK', 'Species', 'YlGn', 1.1)

"""

sorted(pathogensIndexByGenus)
sorted(pathogensIndexBySpecies)

for item in pathogensIndexByGenus:
    print item

for item in pathogensIndexBySpecies:
    print item


#Drawing on Screen
print 'Drawing on screen...'
pl.show()









    
    
