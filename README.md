# PK_goldinglab
A collection of MATLAB scripts, with a small side of documentation, made by Patrick Koch for use in bacterial fluorescent image analysis. 

Note:
These scripts were initially made specifically to be incorperated into the typical workflow in Golding Lab (dated 01/21/2022). Much of the documentation assumes that the reader has prior knowledge of the typical algorithms used in Golding Lab, such as schnitzcell and spatzcell, and their respective and accompanying codes.

Contents of this repository:
	ClusterRecognition
		- This script allows for the identification of closed fluorescent objects and outputs the results in a data format similar to spatzcell.
	
	ClusterAndSpotRecognition
		- This script groups together the spatial information across the spatzcell and cluster recognition data structures. It also includes fake spatzcell objects with randomized spatial information.

	DisplayingSpatzcellImages
		- This collection of scripts, whose original script came from Tianyou Yao, are slight modifications made to better project microscope images across the z-axis into a single XY image.

	GenomeCalculator
		- This MATLAB function takes in E. coli cell growth parameters to output the amount of a specific gene and amount of genome is in the cell during its cell cycle.

	ExaminingImages
		- This script displays the CV, variance, and total intensity of light in each sample’s images and fluorescent channel versus the z-axis. It also plots the z-slice with the highest CV from different fluorescent channels against each other.

	SchnitzcellBadCellMarker
		- This modification of Schnitzcell allows for the user to mark and keep track of specific cells or regions as separate from the background and not a real cell.
