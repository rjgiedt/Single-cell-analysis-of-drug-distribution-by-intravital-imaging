# Single cell analysis of drug distribution by intravital imaging

Installation:

The package contains one code directory and one example dataset directory.

The code is designed to run in Matlab and requires the following Matlab toolboxes:
Image Processing
Statistics
Optimization

The code has been tested Matlab 2017a under Mac OSX 10.7.5. 

To run the code, all files included in the code directory should be placed in the same file location as the video files to be analyzed.  

Running the Pharmacokinetic Data Analysis Program:

The program is designed to be run in the Matlab environment for users with limited knowledge of Matlab.  To run the initial segmentation program open the file “Drug_Analysis.m”.  Clicking the “Run” button in Matlab’s computing environment should then take you through a series of prompts, allowing you to upload your videos, select cell segmentation conditions, and ultimately producing (1) Pharmacokinetic data plots, (2) An output video of cell segmentation selections over time, and (3) an mfile containing a structure for further tracking or linking of your particles for single cell analysis.  Of special importance, the program is formatted to take uncompressed .avi files.  Other file types may result in error generation.    

Single Cell Pharmacokinetic Data Tracking:

As described in the text, users can format the output of the program relatively simply for use in any linking program they are familiar with.  To facilitate single cell tracking, the program has already been formatted to output a structure, “MovieInfo”, that can be uploaded with the linking program “UTRACK”, which is freely available at:

http://lccb.hms.harvard.edu/software.html

To facilitate the combination of the two programs, a script, “Single_Cell_Pharma”, has been provided which can be combined and run with UTRACK, after “MovieInfo” has been generated, allowing for the generation of single cell pharmacokinetic plots with a relatively streamlined process.  For questions regarding single cell tracking, the users are directed to the UTRACK readme.  

Example:
 
One example of drug uptake in an intravital video taken in our lab, is provided to orient users with the algorithm.

To begin analysis, place all files, including codes and video files, into a single directory.  After opening the file, “Drug_Analysis.m” (Green “Go” button in the Matlab environment) you should be prompted to upload two videos; first, a video identifying cell borders or nuclei and second a video of the drug distribution.  After uploading these files, you will be prompted to provide an output file name.

After providing an output video file name, you will be prompted for values to begin segmenting cells.  For this particular example, “good” values are provided along with the description of each parameter.

  
gamma: Is a value controlling the image contrast as is typically seen in imaging applications.  For the example a value of 2 works well.

videoframenumber: Value controls which frame of the video the user would like to analyze.  For the example any value is acceptable, but 1 works well.

iteration: The number of iterations Ray’s method runs through for analysis.  Higher iteration numbers will yield better segmentation at the cost of additional computing time.  For the example 1000 works well.

Power: A measure of the precision of the thresholding method.  For the example, 0.9 works well.

Disk size: A method of removing image speckling, disk size should typically be set to a few pixels.  For the example 3 works well.

Object Minimum: Describes the minimum size of anything that could be an object of interest in pixels.  For the example, 15 works well.  

 Following the selection of pixels, several production images will be produced to consider.  You will then be queried if you would like to try new values (input Y) or not.  If you decline, analysis of the whole movie will proceed.  

The final query in the program asks the user to define the frame rate that they would like to produce an output video at.  Upon completion of the program, several figures, including the average and standard deviation of the drug concentration, minimum and maximum curves, and the minimum drug concentration (if a value is input) number of cells.  

Tracking Cells

To track individual cells, the user is again directed to the UTRACK manual.  The code provided takes advantage of this method to produce plots of all segments analyzed, but analysis of this data is likely to produce more application specific data.

The published manuscript can be found at:
https://www.ncbi.nlm.nih.gov/pubmed/23593370



