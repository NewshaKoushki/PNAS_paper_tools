% TO copy BF images as BigFB images; i.e. ch01 as ch02
% i.e to create a new set of images to fulfil the requirement for a third
% channel 
fdate = 'rawdata';
folds  = {'Cell178';'Cell179';'Cell180';'Cell181' }; % cell folder to analyze

for fi1 = 1: length(folds)      
        cellfolder	= folds{fi1} 
  
       folder1 = strcat('Mark_and_Find_001');
       folder2 = strcat('Mark_and_Find_002');
       
       eval(['cd C:\Singlecell\' fdate '\' cellfolder '\' folder1]); 
       
       d2 = dir; 
       
       for i = 1 : length(d2),
    if  ~isempty(findstr(d2(i).name,'ch01')),	% name must  contain ch01
        copyfile(d2(i).name, strcat(d2(i).name(1:23),num2str(2),d2(i).name(25:28)));

      end;    
    end;
    
     eval(['cd C:\Singlecell\' fdate '\' cellfolder '\' folder2]); 
       
       d2tryp = dir; 
       
       for j = 1 : length(d2tryp),
    if  ~isempty(findstr(d2tryp(j).name,'ch01')),	% name must  contain ch01
        copyfile(d2tryp(j).name, strcat(d2tryp(j).name(1:23),num2str(2),d2tryp(j).name(25:28)));
    end;    
    end;
    
end;