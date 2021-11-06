# SATA

Scripts in this folder is used for the data analysis of LSTM-SOM results. 

## addclicinfo.py

Used to add clinical data in the LSTM-SOM results.

The cluster results were added in`mut_info_XX.csv`

in `addclicinfo.py`, Fisrt the '.csv' files of LSTM-SOM results will be combined, and then devided according to their ICD O3 coding, cancers in some location were merged according to cancer treatment guidelines and clinical practice. Reference: cfg.CFG.clicfeat_dict

## replaceclass.py

Replace the LSTM-SOM class results (i.e. '1210') by new class label (i.e. 'Class_1').



**Before analysis of clinical features, another file of patient information shold be prepared, as the formation of `./SAMPLEDATA/analysis_patient_20200617/thymus.csv `. The example of mutation data is also shown in this folder (./SAMPLEDATA/analysis_data_20200617/thymus.csv)**



## AMBICDclac.py

Used to analyze relationship between MBs classification (datapoint) and clinical features (ICD coding). 

The clinical feature need analysis is assigned in variable 'calcitem' (in __main__ process).

## AMBclicnicalfeaturesata.py

Used to analyze relationship between MBs proportion and clinical features.

## AMBpileclincalfeature.py

Used to analyze relationship between relation between MB classification (datapoint). and clinical features (age, weight, ajcc stage, T stage, N stage, and M stage), the clinical feature need analysis is assigned in variable 'drawfile', and output figure with name 'fig_name' (in __main__ process).

## AMBscattersurvivalforgene.py

Used to draw scatter picture of pairwise survival analysis. Different is reflected by color.

## AMBscartteraspatient.py 

Used to draw scatter picture for MB count and frequency in patients of different cancer (site, pathology).

## Clusterbowmethod.py

Used to decide how many classes is most suitable for K-means cluster. K-means is performed according to the proportion of MBs. A 'elbow' method is used to find out the number of classes. 

## Clusterclinicalfeature.py 

Used to analyze the clinical features of patients in different classes cluster by K-means. The clinical features need analysis are in variable 'CLAC_DICT', if the value is 'avr', it means a continuous variable; if the value is 'chi', it means a discrete variable. 

## Clustersurvivalanalysis.py

Survival analysis for patients according to classes cluster by K-means.

## Drawbasepile.py 

Draw histogram for base proportion beside mutation site.

## gene_sata.py

Draw heatmap for genes with high frequency in different classes of MBs.

