November 18th, 2016
Cell Tracking Code

Updated by Sam Bruce, Ajay Gopinath, Kendra Nyberg, Mike Scott 
Original Code by Bino Varghese, updated by David Hoelzle

Code is designed to track cells in the cell deformer device and output their transit times and areas.  It is designed to have minimal user input so reliable results are obtained.

To use this code:
	- Ensure video files are in .avi format, and named correctly.
	For a 5 micron device at 800 fps: 
		'dev5x10_800fps[anything you want].avi'
	- Run MainCode.m
	- Select user preferences when prompted
		- GUI will pop up to input user preferences
		- Select the type of particle to track (e.g. cell, gel or oil particle)
		- Select whether or not to check for correct mask fitting
			- Yes: data is not processed and a mask overlay is displayed
			- No: data is processed and assumes correct mask fitting
		- Select if there are manual detection filters.
	- Select videos when prompted
		- GUI will pop up to select videos from a folder
		- Select as many video as you want from one folder
		- GUI will pop up to select videos from other folders
	- Code will automatically detect constriction size and framerate
	- A progress bar will tell you how much longer the code will run

Output:
	- An excel sheet of data for the first constriction transit time data
	- Each row represents a single cell/particle
	- Each column contains data for the following:
		1.  Effective Diameter
		2.  C1 transit		3.  Entry Frame		4.  Exit Frame		5.  Total Occlusion at entry		6.  Total Occlusion at exit		7.  Average Occlusion during transit		8.  Maximum Occlusion during transit		9.  Lane
		10. Parent directory
		11. Daughter directory
		

This version of the code has many updates from the original code, many of which are listed below.
	- Constriction size and framerate are detected from the filename
	- Multiple select is enabled
	- Individual processed frames are not saved to the hard drive
	- Digital image filtering was improved
	- Multiple average backgrounds are now calculated
	- Searching algorithm was redone, and no longer searches
	- Preprocessing of 50 frames was eliminated
	- Calls to regionprops were drastically reduced for speed
	- Cell area is calculated before entry into the lane and at 				constrictions
	- Cells are separated into lanes for analysis
	- Masks are used to define the channels
	- Autocorrelation with the mask is used to define constriction 				regions
	- Videos are no longer cropped (didn't significantly alter speed)
	- 'Paired cells' (those where 2 cells are in the same lane 				concurrently) are separated from cells who pass through an 				otherwise empty lane.
	- Progressbar labels are now persistant
	- Cells that touch multiple lines are now counted at the lowest 			(highest numbered) line they touch.  Solves timing issues 				caused by the line used to calculate area.
	- Detection of gel and oil particles is enabled
	- Inputs of manual detection filters is included as an option
	- Transit times for all subsequent constrictions (C2-C6) are removed from data 		output.
	- Lane occlusion is added to analysis to keep track of occupied neighboring 		lanes.

