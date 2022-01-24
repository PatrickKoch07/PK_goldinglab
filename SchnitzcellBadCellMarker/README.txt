General idea:
The point of this edit in schnitzcell is to allow for the user to mark specific cells as special.

The edit does this by simply making that desired cell's cell number negative. By designating some cells, or regions if the user likes, as distinguishable from other cells, we can do analysis in spatzcell and in FluorMeasure without including possible experimental defects. For example, we can mark overlapping cells to not be included as cells for analysis and not relegate them to the background. Another example can be cells which are divided by the image boundary. These cells should not be included for analysis nor should they be a part of the background.


How to run the code:
SchnitzcellBadCellMarker_output.png
To use this code, run schnitzcell with the appropriately changed code attached in the previous page. In the console there is an instruction printed which tells the user which button to press to mark the cell. In the image above we can see this console display and we see that after applied to the cells on the boundary, they are marked.

It's important that in any follow up code that the user runs that they change the behavior of the script to analyze positive cell numbers only.

Example:
	LcFull = uint16(LcFull);
	T_LcFull = uint16(zeros(size(LcFull)));

These two lines in spatzcell will get LcFull and convert them to be unsigned integers (turning them into zeros). This should effectively make the new LcFull a collection of only the cells we want analyzed.
