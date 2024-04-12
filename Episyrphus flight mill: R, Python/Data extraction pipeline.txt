The data extraction pipeline - order of files

Raw data (inputed manually):
Data S2/Fly information.csv - Information for each fly (condition, sex, size etc)
Data S2/Morphometric information.xlsx - Raw morphometric data

Raw data (created by Code S3/1_Flight_mill_running_code.py):
Data S2/21 Autumn/Mills/[mmdd] - Raw data for flight mill experiments on autumn migratory population
Data S2/22 Summer/Mills/[mmdd] - Raw data for flight mill experiments on summer non-migratory population

Processed data (created by Code S3/2.1_Data_extraction.py):
Data S2/[Experiment season]/All flights.csv - Flight metrics on individual flight events
Data S2/[Experiment season]/Flight_stats.csv - Flight metrics per fly
Data S2/[Experiment season]/Transects.csv - Time transect (requires aggregation before analysis)
Data S2/[Experiment season]/All flights.csv - Information on individual flight events

Aggregated processed data (created by Code S3/2.2_Data_concatenation,R):
Data S2/RD All flights.csv - Flight metrics on individual flight events
Data S2/RD Flight_stats.csv - Flight metrics per fly
Data S2/RD Transects.csv - Time transect
Data S2/RD All flights.csv - Information on individual flight events

Statistical analysis and graphs (created by Code S3/3.1_Mill_Analysis.R)
Morphometric analysis and graphs (created by Code S3/3.2_Morphometrics.R)

##############
Morphometric information

Size:
Manuscript levels = small, moderate, large.
Data & analysis code levels = small, average, large
