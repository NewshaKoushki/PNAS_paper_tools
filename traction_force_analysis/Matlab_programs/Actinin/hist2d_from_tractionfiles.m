% This code uses the traction .dat files to plot the evolution of the
% distribution of the traction values in time, as a 2D histogram. 
% Matrices necessary to plot the 2d histrogram are save in histo.mat in the
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
    
    % Create a big time matrix
    if exist('time_SumUp','var') == 0, % initialise the variable at the 1st step of the loop
    time_SumUp = zeros(length(x),time_points);  
    end
    % Fill in the SumUp matrix with time amplitudes
    time_SumUp(:,t)=t;
    
   
    
end %loop over time points 
   
% PART FOR REGULAR 2d histogram
T_SumUp_v = reshape(T_SumUp,[],1);
idx_zeros = find (T_SumUp_v <5);
T_SumUp_v(idx_zeros) = NaN;

time_SumUp_v = reshape(time_SumUp,[],1);

data=cat(2,time_SumUp_v,T_SumUp_v);
bin_time = (2.5:5:142.5);
bin_force = (100:200:6900);
n = hist3(data, {bin_time,bin_force});

nH1 = n';
nH1(size(n,2)+1,size(n,1)+1) = 0;
xb = linspace(0,150,size(n,1)+1);
yb = linspace(0,7200,size(n,2)+1);
figure,
h = surf(xb,yb,-log((nH1)/sum(sum(nH1))),'EdgeColor', 'none'); 
colorbar
axis([0 145 0 7000])
colormap (hot)
view(2)
box off 
grid off
caxis ([4 11]);

save('histo.mat','xb','yb','nH1')

clear time_SumUp time_SumUp_v T_SumUp T_SumUp_v


 end % end loop over folders
 
 
 
 
 
 
