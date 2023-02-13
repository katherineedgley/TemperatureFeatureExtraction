# TemperatureFeatureExtraction

Source code for extracting temperature features as reported in the study "New Insights into Stroke from Continuous Passively Collected Temperature and Sleep Data Using Wrist-Worn Wearables", which can be found here: https://www.mdpi.com/1424-8220/23/3/1069

****************************************


## **Instructions for use**   

**Extraction from single file** (_single_file_extraction.R_): 

Features are extracted for a single ".bin" file ("filename") using the function _get_temp_features_. Requires output from GGIR (here using version 1.9-2, configuration as reported in article) to use sleeping times, stored in "output_dir". Later versions of GGIR are not compatible - we are working on updating this code to keep it compatible with the most recent versions of GGIR as well. 

Indicative call: get_temp_features("patient_0.bin", "/Desktop/GGIR_output")

##


**Extraction from multiple files** (_multiple_file_extraction.R_):

Multiple files can be processed at once using the _process_file_ls_ function. As input it requires a vector of ".bin" filenames (output from GENEActiv watches), and the "output_dir" where GGIR (version 1.9-2) output is stored. 

Indicative call: dt <- process_file_ls(c("patient_0.bin", "patient_1.bin"), "/Desktop/GGIR_output")


****************************************

## **If you use this code please cite:**

Edgley, K.; Chun, H.-Y.Y.; Whiteley, W.N.; Tsanas, A. New Insights into Stroke from Continuous Passively Collected Temperature and Sleep Data Using Wrist-Worn Wearables. Sensors 2023, 23, 1069. https://doi.org/10.3390/s23031069 
