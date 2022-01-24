How to run the code (going line by line down the code for things to change):

	- Change the function whose return value is stored in 'ip'. This should be at the start and called "exp_YYYYMMDD_InitializeExp;"
	- Change all channel information such as color, loading information, and data folders
	- Change sample names
	- Set random seed as desired (for when we select random spot locations)
	- Set desired default minimum distance if there are no objects to which to find a minimum distance.
	- Set desired radii for checking cluster fluorescence and radii for checking channel 2 spots.
	- Again make sure segmentation folder is correct for cell masks
	- Again make sure that the datafolder for peakdatas are correct
	- Again make sure cluster mask folder and cluster fluorescence image folder are correct
	- Within "InitializeClusterRecognitionParameters_c234" make sure the parameters are correct for this experiment (only used to load in the cluster fluor image stack)
	- There is one commented out paragraph of code which allows you to visualize the overlap between the channel 1 spot and cluster masks. If this is uncommented, it would be wise to put a breakpoint here so many figures aren't left open.
	- Check at the end of the code if this is where you want the data to be saved

What to expect as outputs:

As output a data set, distancelist_new, will be made and saved:
	ClusterAndSpotRecognition_output1.png
Ex. of one row of data

There should be 4 rows of data:
	Channel 1 real spots (labeled as Cy3 above)
	Channel 1 fake spots
	Channel 2 real spots (labeled as Cy5 above)
	Channel 2 fake spots

The columns of the data set will be the samples.

For each element in the data set (so for a particular channel and spot type, and for a particular sample) there should be 21 columns of data as seen above with each row corresponding to a spot.
