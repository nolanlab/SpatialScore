# SpatialScore
Code for calculating the SpatialScore as described in our manuscript: "Immune cell topography predicts response to PD-1 blockade in cutaneous T cell lymphoma"

Introduction
Using the X/Y coordinates for each cell-type, this R application calculates the spatial distances (i.e., the minimal distance between each cell-type and its nearest other cell-types) as well as ratio distances between 3 cell-types.

In our manuscript, we were specifically interested in the distance relationships between 3 cell-types: effector T cells (CT1), tumor cells (CT2) and Tregs (CT3). We used this R application to calculate the ratio of the minimal distances between CT1—CT2 (right distance) versus CT1—CT3 (left distance). This distance ratio represents the SpatialScore.

This application also assess whether these distance ratios are significantly different from those of a random sample. For the number of CT1 cells in each tissue region, we randomly selected the same number of non-CT1 (nCT1) cells. For each of these nCT1 cells, we calculated the ratio of the minimal distances (nCT1—CT2 / nCT1—CT3) and determined the mean of this sample. We repeated this random sampling 100 times, and the average of all the means was reported. Distribution of the random values was assessed by the quant output variable, which indicates how many of the random means are smaller than the measured means. For instance, a quant of 97 indicates that 97% of the random means are smaller than the measured means. Thus, quant values closer to 100 or 0 indicate that the measured means are not random.

Installation: this is an R script, so it just needs downloading. The script will install all required packages on the first run.

Input: data table (txt) containing several required pieces of information for each cell:

  - x, y coordinates

  - tissue section id

  - cluster assignment

  - group assignment

Output: large number of intermediate files (output from Delauney algorithm, plus processing files that can be ignored):

  - heatmaps for niches per cell type

  - heatmaps for cell types vs cell types per groups using likelihood ratios (llr)
  
  - if applicable: comparison of llr between groups
  
  - density plots or histograms of minimal distances for each cell type pair – one density/histogram for each group.

  - data table of minimal distances to each cell type per each cell (2 files: in pixels and um, here assuming a 0.37744 conversion factor based on our objective)

  - data table of ave minimal distances and their statistics, per spot

  - data table of minimal distance ratios and their statistics, per spot
