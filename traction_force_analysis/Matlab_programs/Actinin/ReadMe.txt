Notes to run the codes.

1/ you need to copy the 'SingleCell' folder on the C: drive of your PC so that the paths used in the codes work. 

2/ In the 'Matlab_programs' folder, you will find

area_centroid.m
% This code uses the files containing the cells boundaries (displacement folder)
% to create 3 .mat files containing:
% 1/ the evolution of cell area 
% 2/ the trajectory of the centroid
% 3/ the msd of the centroid  

hist2d_from_tractionfiles.m
% This code uses the traction .dat files to plot the evolution of the
% distribution of the traction values in time, as a 2D histogram. 
% Matrices necessary to plot the 2d histrogram are saved in histo.mat in the
% corresponding cell folder of the traction folder.

hotspots_from_tractionfiles.m
% This code uses the traction .dat files to look for the hospots defined as 
% traction>0.5*traction_max, adds up the hotspots over time to get an idea 
% of their persistency 
% Matrices necessary to plot the results in a histrogram are saved in histo_cum.mat in the
% corresponding cell folder of the traction folder.

importfile.m is an 'homemade' function used in some of those codes. 

You can not use these codes if you don't have the files containing the cell boundaries or the traction files (.dat)
I didn't upload them because they are too big but I thought sending you these codes would be useful. 

3/ The outputs I got with these codes are saved as .mat files in the 'displacement' (msd and cell area) and the 'traction' folders (traction 2d histograms and hotspots info). 
I sorted them in groups (manually). 

4/ In the 'Matlab_programs' folder, you will also find 

plot_msd_area.m
% This code loads the .mat files that gather the msd and the area info for each cell of
% the different groups, and plot them.

plot_final_hist2d.m
% This code loads the .mat files that gather the 2d histogram info for each cell of
% the different groups, and plot one 2d histogram per group.

plot_final_histcum.m
% This code loads the .mat files that gather the histogram info to represent 
% the hotspots persistency for each cell of the different groups, and plot 
% 2 histograms to compare the groups at different stiffnesses. 

You can run these codes to obtain the figures I copied in the ppt (except the hotspots maps), 
And you can modify them if you want to play with the representations. 

