# SpatialScore
Code for calculating the SpatialScore as described in our manuscript: "Immune cell topography predicts response to PD-1 blockade in cutaneous T cell lymphoma"


**Introduction**

Using the X/Y coordinates for each cell-type, this R application calculates the spatial distances (i.e., the minimal distance between each cell-type and its nearest other cell-types) as well as ratio distances between 3 cell-types.

In our manuscript, we were specifically interested in the distance relationships between 3 cell-types: effector T cells (CT1), tumor cells (CT2) and Tregs (CT3). We used this R application (spatial_analytics.R) to calculate the ratio of the minimal distances between CT1—CT2 (right distance) versus CT1—CT3 (left distance). This distance ratio represents the SpatialScore.

![SpatialScore](https://user-images.githubusercontent.com/37353112/134821486-15bd8d31-a134-4143-8851-b96a85ff5292.PNG)

This application also assess whether these distance ratios are significantly different from those of a random sample. For the number of CT1 cells in each tissue region, we randomly selected the same number of non-CT1 (nCT1) cells. For each of these nCT1 cells, we calculated the ratio of the minimal distances (nCT1—CT2 / nCT1—CT3) and determined the mean of this sample. We repeated this random sampling 100 times, and the average of all the means was reported. Distribution of the random values was assessed by the quant output variable, which indicates how many of the random means are smaller than the measured means. For instance, a quant of 97 indicates that 97% of the random means are smaller than the measured means. Thus, quant values closer to 100 or 0 indicate that the measured means are not random.


**Installation**

This is an R script, so it just needs downloading. On the first run, the script will install (or prompt you to manually install) any required packages. 


**Input and Output Data**

Input: a data table (.cvs) containing several required pieces of information for each cell:

  - Cluster Column Name: the column name that contains the cell-type (e.g., Tregs, tumor, etc.)
  - Patient/Treatment Group Coulmn Name (e.g., column name that contains responder=1 vs non-responder=2)
  - Column Name for Regions (e.g., column name that contains tissue microarray spot number or different tissue regions)
  - Column Name for Patients (e.g., column name that contains the patient ID)
  - X Coordinate Column: column that cotains the X postion for each cell
  - Y Coordinate Column: column that cotains the Y postion for each cell
  - Number of niches: not relevant for current study (default value = 10)
  - Minimum cell count total: default value = 100
  - Conversion factor (um/pixel): dependent on the objective used for imaging (0.37744 for this study) 

Output: various data tables (.txt) and numerous intermediate/processing files that can be ignored

  - Files used for heatmaps and niches per cell-type
  - Heatmaps for cell types vs cell types per groups using likelihood ratios (llr) and if applicable, comparison of llr between groups
  - Density plots or histograms of minimal distances for each cell type pair (one density/histogram for each group)
  - Data tables of minimal distances to each cell type per each cell (2 files: in pixels and um)
  - Data tables of minimal distances and their statistics, per spot and per group
  - Data tables of minimal distance ratios and their statistics, per spot and per group
  - Data tables comparing measured vs random distance ratios using Method 1 and Method 2. Method 1 was used in this study; z-score was calculated as follows: z-score = (ratio_mean - mean_rand) / sd_rand


**Test Data**

Here we provide an sample input data table (SpatialScore_test_data.csv) and the resulting output data tables for reference. 


**Script creation**

This script was created by Janos Demeter, PhD (janos.demeter@stanford.edu)  
