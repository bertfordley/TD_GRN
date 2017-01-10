Author: Rob Moseley
Email: rmosele4@vols.utk.edu

This was created to create a gene regulatory network from time course data.
This is written specifically for my data (e.g. hard coded areas) but would like to make
it available for any time course data.

Input Arguments:
	:input_file         data file (see Input Data Format section below for appropriate formatting)
	:s                  max value for linearly spaced vector (for cubic spline interoplation)
	:r                  max time delay
	:thres              threshold for skipping genes in time sample alignment (for cxcorr values)
Input Data Format:
    - 3 bioreps per 12 time point samples
    - rows = genes
    - columns = time points
Steps:
    - computes averages for each time point sample
    - converts data to 1-hour intervals
    - converts data to z-scores
    - runs cross correlation to find time delays
    - align time samples and perform spearman rank then network deconvolution for each gene
Outputs:
    - Averaged data
    - Z-Scores
    - Direct dependencies matrix in json format
Usage:
    - python CcorrND.py <input_file> <s> <r>

'''