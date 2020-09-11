# photometry

-When you download these files, put these files in your MATLAB path:
getStates
plotStates
group_z_project_vector
baseline1
segmentFunction
getStateAvgs

- Put these files in the folder with your data:
fancyPlot_portable
segmentData
group_transitions
quantifyStagesMay20
multiMouseAvgs

- the yew files are example data sets that work with these scripts

1.	Check for noise in the recording:

To remove noise, run these commands (will introduce 'nan' values into your data):

- Mouse = 'yew1501.mat';
- load(Mouse);
- plot(Calcium.values);
*** Then look at the plot and see where the noise is, exclude all data points that are less than (or greater than, depending on the noise) a value that includes all the real data
- Calcium.values(Calcium.values<0.67) = nan;

Then use ‘Ctrl click’ to select Calcium and instant_state from the Workspace then right click to save these variables as a new ‘.mat’ file to analyze your data with. 


2.	Run fancyPlot_portable to see individual data:
-	Needs scripts ‘getStates’ and ‘plotStates’
-	You might need to fit a modified curve to the data. Use the curve fitting app to do this.
-	Remember to save your fitted model as ‘coeff_mouseName.mat’ (where mouseName is the name of your datafile)

3.	Run segmentData:
-	Needs scripts ‘group_z_project_vector’, ‘baseline2’ and to ‘segmentFunction’ work
-	Select all your options for the segmentData script
-	Remember that this script will execute on ALL FILES in the current working directory that begin with whatever you put in the variable ‘names’.
-	Run this script for each different state transition (i.e. I usually run it 4 times, changing settings.stage and settings.stage2 so I get:
1.	NREM-REM
2.	REM-WAKE
3.	WAKE-NREM
4.	NREM-WAKE
  
4.	Run group_transitions script to average over mice

-	sleep wake states (from vs. to) in opposite order to that written in the segmentData script
-	need to change cmin, cmax, ymin, ymax so all graphs are plotted on same (best) scale
-	the pre_state_time and post_state_time MUST be set to the same values as specified in the segmentData script (the groupTransitions script uses the segmented data from the segmentData script)
-	time_value_pre and time_value_post will average the signal for a specified time before and after the transition point. I usually set this to 10 seconds each side and run my statistics on the results to determine whether there is a significant difference in activity before and after a transition. The result is calculated and stored in the variables 'meanMinPre' and 'meanMinPost'. You can select these variables by double clicking on the variables name Workspace window and then copy and paste the values into GraphPad to run t-tests and things. 

5.	Average net GCaMP6s activity by sleep-wake state
a.	Run quantifyStagesMay20 for each mouse  (this script requires ‘getStateAvgs’)
b.	Run multiMouseAvgs
i.	Select the ‘averages’ folder from the working directory
ii.	Copy wholeDataSet to your stats package. The rows are each mouse, the columns are the average wake, NREM and REM value for each mouse



To export any graphs as an emf 
1.	Click on the figure you want to save (do not click any other window before doing the next step)
2.	Run this command in the command window:
3.	set(gcf,'renderer','painters')
4.	Go back to the figure and click the save icon on the figure. In the 'save as type' option, select 'Enhanced metafile (*.emf)'
