How to run the code (going line by line down the code for things to change):

	- Change the function whose return value is stored in 'ip'. This should be at the start and called "exp_YYYYMMDD_InitializeExp;"
	- Tell the code which fluorescent channels are which channel index.
	- Make sure 'norm_over_all' is set to 0. This way each cell is normalized over itself rather than with respect over the whole fluorescent image.
	- Set the 'combo' of 'high_threshold' and 'low_threshold' to be run. The reason why a histogram is produced for each image is so that you can go back and decide for yourself a good value of high and low threshold to place here. This also makes it easy to test a few thresholds all at once.
	- Make sure the 'spot_folder' is how you wish for it to be named/organized.
	- Check the function file "InitializeClusterRecognitionParameters_c234" to make sure it points to the correct file locations for the microscope and segmentation images.
	- Make sure the disk radii for the 'imerode' and 'imdilate' functions are what you want them to be.

What to expect as outputs:
	
	As mentioned previously, for each image we will be returned a histogram of pixel intensities with a three gaussian curve fit. This three gaussian curve fit is supposed to approximate a three pixel intensity population within the cell: Background, low threshold area, and high threshold area. 
	
	For each image we will also be returned a reconstruction of the normalized fluorescent image, the threshold objects without filters, and the final filtered objects. More detail into the displaying algorithm can be found in the 'Displaying Images' tab.
	
	Finally, we are also returned data regarding the threshold objects:
		- A matrix of enumerated low threshold objects is stored for each image (I have elected to REMOVE the high threshold areas, this can be changed by commenting out the following line:
			§ " Lcluster_low = Lcluster_low .* (Lcluster_low & ~smooth_T_high_image); "
		- A matrix of enumerated high threshold objects is stored for each image.
		- Below is what is returned for the low threshold objects and how they are extracted from the images:
			§ ClusterRecognition_output1.png
		- Below is what is returned for the high threshold objects and how they are extracted from the images:
			§ ClusterRecognition_output2.png
			§ Note that the only difference is the addition of the final column labeling the high cluster.
