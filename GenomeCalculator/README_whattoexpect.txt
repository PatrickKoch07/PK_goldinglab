How to run the code:

To run this code just follow the inputs of the function as explained in the "Genome calculator" tab along with selecting your gene's name in gene_name and telling the code if you wish to output figures as well with show_fig. The C and D period can be estimated, so they don't have to be given. However note that this is an estimation for K-12 MG1655 from Michelsen, et al., "Precise determinations of C and D periods by flow cytometry in Escherichia coli K-12 and B/r". Microbiology (2003):
	
	GenomeCalculator_C&D_prediction.png



What to expect as outputs:

For outputs, there will be 2 vectors and a scalar returned by the function. geno will hold genome copy number vs. the requested time points, gene will hold gene copy number vs. the requested time points, and gene_time will give the time from birth that the gene replicates.

If show_fig is true, then the following figures will be outputted:
	GenomeCalculator_output1.png