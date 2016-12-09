Author: Rob Moseley
Email: rmosele4@vols.utk.edu

This was created to create a gene regulatory network from time course data.
This is written specifically for my data (e.g. hard coded areas) but would like to make
it available for any time course data.


Input Data Format:
    - 3 bioreps per 12 time point samples
    - rows = genes
    - columns = time points
Steps:
    - computes averages for each time point sample
    - converts data to 1-hour intervals
    - converts data to z-scores
    - runs cross correlation and network deconvolution for each gene
    - ranks gene-gene interactions
Outputs:
    - Averaged data
    - Z-Scores
    - Direct dependencies matrix
    - Time Delay matrix
    - Ranked List of gene-gene interactions
Usage:
    - python CcorrND.py <input_data> <spline_number> <max_time_delay> <maxcount>
