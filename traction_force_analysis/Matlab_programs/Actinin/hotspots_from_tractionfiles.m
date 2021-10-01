% This code uses the traction .dat files to look for the hospots defined as 
% traction>0.5*traction_max, adds up the hotspots over time to get an idea 
% of their persistency 
% Matrices necessary to plot the results in a histrogram are saved in histo_cum.mat in the
% corresponding cell folder of the traction folder.


clear all

%Mutants 4kPa
folds  = {'Cell3101','Cell3102','Cell3103','Cell102','Cell3110'}; % cell folder to analyze 
time_points = 51;

%Mutants 26kPa
%folds  = {'Cell103','Cell3113','Cell111','Cell3117','Cell3118','Cell3120'}; % cell folder to analyze 
%time_points = 145;

%Ctrl 26kPa
%folds  = {'Cell104','Cell3123','Cell3125','Cell112','Cell3129'}; % cell folder to analyze 
%time_points = 145;

%Ctrl 4kPa
%folds  = {'Cell105','Cell3132','Cell3134','Cell3136','Cell3137'}; % cell folder to analyze 
%time_points = 51;

 for fi1 = 1: length(folds)      
        cellfolder	= folds{fi1}
eval(['cd C:\Singlecell\traction\' cellfolder]); 

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% FIRST, SECOND - counters needed to get the right number of time points based on the
% number of files in the folder

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
 
posn = 1; % Position can be put in a loop if needed

%=========================================================================
% Analyse Hotspot persistency over time

for t=1:time_points, 
    
    % Getting the name of the data file (time point)
            
    if  posn <= 9,
                  if t <= 10,
                  name = strcat('Pos00',num2str(posn),'_S001_t00',num2str(t-1),'_1.dat');
                  else
                      if t <= 100,
                      name = strcat('Pos00',num2str(posn),'_S001_t0',num2str(t-1),'_1.dat');
                      else
                      name = strcat('Pos00',num2str(posn),'_S001_t',num2str(t-1),'_1.dat');
                      end;
                  end;
    else
                  if time <= 10,
                  name = strcat('Pos0',num2str(posn),'_S001_t00',num2str(time-1),'_1.dat')
                  else
                      if time <= 100,
                      name = strcat('Pos00',num2str(posn),'_S001_t0',num2str(time-1),'_1.dat')
                      else
                          name = strcat('Pos00',num2str(posn),'_S001_t',num2str(time-1),'_1.dat')
                      end;
                  end;
    end
        
        
        
    % Import traction data file
    
    [x ,y ,tx, ty] = importfile(name,2,Inf);
    clear fi1; %clear the file identifier after importing a file to avoid mixing up during next use

    
    % Create a big T matrix that gathers all columns 
    if exist('T_SumUp','var') == 0, % initialise the variable at the 1st step of the loop
    T_SumUp = zeros(length(x),time_points);  
    end
    
    % Fill in the SumUp matrix with traction amplitudes
    T_SumUp(:,t) = sqrt(tx.^2+ty.^2);
    
    % Create a big H matrix that gathers info about Hotspots 
    if exist('H_SumUp','var') == 0, % initialise the variable at the 1st step of the loop
    H_SumUp = zeros(length(x),time_points);  
    end
    
    % Fill in the SumUp matrix with Hotspots (1) or not (0) 
    Tmax = max(T_SumUp(:,t)); 
    idx = find(T_SumUp(:,t)>0.5*Tmax);
    H_SumUp(idx,t) = 1;
   
end %loop over time points 

% PART FOR REGULAR histogram
H_Cum_v = sum(H_SumUp,2);
idx_zeros = find (H_Cum_v == 0);
H_Cum_v(idx_zeros) = NaN;
bin_time = (2.5:5:142.5);
hist_cum = hist(H_Cum_v,bin_time);


save('histo_cum.mat','bin_time','hist_cum')

clear T_SumUp H_SumUp
 end % end loop over folders

