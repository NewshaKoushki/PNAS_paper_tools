% This code uses the files containing the cells boundaries (displacement folder)
% to create 3 .mat files containing:
% 1/ the evolution of cell area 
% 2/ the trajectory of the centroid
% 3/ the msd of the centroid 

% INPUT: folder name (in the displacement folder) and pixel size 
% OUTPUT: .mat files saved in the cell folder of the displacement
% folder 

close all;
clear all;
     
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%  INPUT: PIXELSIZE, YOUNG'S MODULUS, POISSON'S RATIO


folds         = {'Cell3102'};
ncells       = 1;
pixelsize     = 0.2522;

for fi1 = 1: length(folds)      
        cellfolder	= folds{fi1}
        
        
eval(['cd C:\Singlecell\displacement\' cellfolder]); pwd; 
        
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% FIRST, SECOND
   d2 = dir; si = 1;
   d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
   d2b(2)=d2(2);
   for d2i = 3 : length(d2),
      if ~isempty(findstr(d2(d2i).name,'.dat')) &  ...		% search for an excel file
          ~isempty(findstr(d2(d2i).name,'Pos')), 			% name has to contain 'trypsindisp'
      		d2b(si+2)	=d2(d2i);
		    d2b(si+2).name 	= d2(d2i).name;
            si = si + 1;
      end;
   end;

   
   
for kloop = 3: length(d2b)

    eval(['cd C:\Singlecell\displacement\' cellfolder]); pwd; 
  
    savefilename = d2b(kloop).name; % needed for saving figures    
    displ = load(d2b(kloop).name);
    digit = findstr(savefilename, '_ch');
    
if exist('area_cell','var') == 0,
    area_cell = zeros(length(d2b)-3,1);   % initialise the variable at the 1st step of the loop
end  

if exist('centroid_cell','var') == 0,
    centroid_cell = zeros(length(d2b)-3,2);   % initialise the variable at the 1st step of the loop
end  
 
%====================================================================================

 for   cellno = 1 : ncells,  % loop over number of cells per image 
        
        
      
      
   filename = strcat(savefilename(1:digit-1),'_',num2str(cellno),'.mat');
   load(filename);
   xrub = (xrub - mean(displ(:,1))) * pixelsize;
   yrub = (yrub - mean(displ(:,2))) * pixelsize;
 
        eval(['cd C:\Singlecell\Matlab_programs']); pwd;
   [centroid_cell(kloop-2,:), area_cell(kloop-2)] = polygonCentroid(xrub,yrub);
        
 end % number of cells per image (here =1)       
 
end % kloop time

%% calculate MSD 
%(taken and adapted from http://stackoverflow.com/questions/7489048/calculating-mean-squared-displacement-msd-with-matlab)

nData = size(centroid_cell,1); %# number of data points
numberOfDeltaT = floor(nData); %# for MSD, dt should be up to 1/4 of number of data points

msd = zeros(numberOfDeltaT,3); %# We'll store [mean, std, n]

%# calculate msd for all deltaT's

for dt = 1:numberOfDeltaT
   deltaCoords = centroid_cell(1+dt:end,:) - centroid_cell(1:end-dt,:);
   squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2

   msd(dt,1) = mean(squaredDisplacement); %# average
   msd(dt,2) = std(squaredDisplacement); %# std
   msd(dt,3) = length(squaredDisplacement); %# n
end

eval(['cd C:\Singlecell\displacement\' cellfolder]); pwd;  

%% save main data
save('area_cell.mat', 'area_cell')
save('centroid_cell.mat', 'centroid_cell')
save('msd.mat', 'msd')

end % cell folder



