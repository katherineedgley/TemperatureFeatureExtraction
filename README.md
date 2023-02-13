# TemperatureFeatureExtraction

Source code for extracting temperature features as reported in the study "New Insights into Stroke from Continuous Passively Collected Temperature and Sleep Data Using Wrist-Worn Wearables", which can be found here: https://www.mdpi.com/1424-8220/23/3/1069

****************************************


**Instructions for use** 

Features are extracted for a single ".bin" file ("filename") using the function _get_temp_features_. Requires output from GGIR (here using version 1.9-2, configuration as reported in article) to use sleeping times, stored in "output_dir". Later versions of GGIR are not compatible - we are working on updating this code to keep it compatible with the most recent versions of GGIR as well. 

_get_temp_features_ can be modified to iterate through multiple ".bin" files as indicated in function. 

****************************************

**If you use this program please cite:**

Edgley, K.; Chun, H.-Y.Y.; Whiteley, W.N.; Tsanas, A. New Insights into Stroke from Continuous Passively Collected Temperature and Sleep Data Using Wrist-Worn Wearables. Sensors 2023, 23, 1069. https://doi.org/10.3390/s23031069 
