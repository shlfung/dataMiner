For any question, ask Stephen Fung (stephen.fung@agr.gc.ca) at Pacific Agri-Food Research Centre in Agassiz, BC

This script is primarily for local computer (where all the experimental results are stored the local computer) with limited memory. A version that will be optimized for the biocluster computer is in development

Instruction for using the dataMiner script on a local desktop:

The script will search and analyze the results of a given organism from the Automated Identification Pipeline (Auto ID) and the Oligo Fishing Pipelines (OFP) of the SAGES project. The output is a csv file that indicates that abundance of the bacteria/fungi (input by the user) in each sample with the relevant meta information (Location, collector, date, climate…etc) of each sample. The csv file can be used for downstream analysis (development in progress)

It requires Python 2.7 with Biopython and Pandas installed 

To execute the script, in command line:

python dataMiner.py "rank of the organism classification" "name of the species" "Directory of the Summary Logs (Auto ID Results) with the file extensions" "SAGES Meta information in CSV"  "OFP Result File" "Yes/No"

Yes = Display ONLY the Samples with the organism in the output
No = Display all the Samples in the output

Example:

python dataMiner.py genus Streptococcus ./rs75_Summary_Log/*.summaylog ./input/454_sample_summary_with_climate.csv ./input/Streptococcus_spp.Streptococcus_spp.txt.withmeta.fna.fasta No

This will produce a csv file that indicates which sample contains the Streptococcus bacteria. The abundance of Streptococcus in each sample will be displayed. The meta information (location, collector, date, climate..etc) will also be displayed in the csv file
