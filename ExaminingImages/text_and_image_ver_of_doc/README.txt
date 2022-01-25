General idea:
To examine our images, we'll look at an images intensity, variance, and CV. There exist three versions of this, one for before segmentation and two for after segmentation. This script simply finds these values from the given image set and plots the results.

	- By looking at an image's intensity vs. z-slice we can get a sense of the amount of bleaching occurring in the image.
	- By looking at an image's variance we can get a sense of the crispness of edges.
	- By looking at an image's CV, we can hopefully gain the benefit of looking at variance without any side effects of bleaching, or changing intensity values between z slices.

	- Note that making these plots over a whole frame can be drastically different from making these plots over each cell in a frame. If there is a drastic difference then it will be known that analysis which uses whole frame focus might need to be changed.

	- From phase 3, by looking at the CV vs. z-slice we can see if our image's best focus, as determined by the highest CV value, lies at least between the two end slices. 

	- Next, by plotting a fluorescent channel of interest's best slice vs. phase 3's best slice we can see how closely the two match. If, as verified by also going back to look at the microscope images, there is no similarity of focus between the two channels, then code where phase is used to focus the fluorescent channel should be changed or not used.

		? Additionally this plot can show if the average fluorescence best slice is consistently focused above or below the phase 3 channel. Using this, after checking microscope images, can help determine if phase 3 can be used to focus the fluorescent channel.


The following plots are examples of what is created:
	ExaiminingImages_output1part1.png
	ExaiminingImages_output1part2.png
	
	ExaiminingImages_output2.png

	ExaiminingImages_output3.png