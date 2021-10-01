% Program to extract specific data entries from an excel file
% Written in order to extract specified displacement values in frames
% closest to magnetic bead position
% To test whether synchronous displacements around the bead change in time
% To check whether there are other no synchronous displacement maxima
% Written by Ramaswamy Krishnan on 11/1/06
%  INPUT: LOAD DISPLACEMENTS

%====================================================================================    
   close all;
   clear all;
   
     folds = {'Cell27','Cell28','Cell29'};
     matrix = [];
     
for fi1 = 1 : length(folds),
    cellfolder	= folds{fi1}
    eval(['Syn_mvmt']); 
    eval(['cd C:\C2C12\Displacement\Syncronized\' cellfolder]);
    pwd;
    savedirectory = strcat('C:\C2C12\traction\Syncronized\',cellfolder);

%============================================================================================
% Identify displacement files
   d2 = dir; si = 1;
   d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
   d2b(2)=d2(2);
   for d2i = 3 : length(d2),
      if ~isempty(findstr(d2(d2i).name,'.xls')) & ...		% search for an excel file
          ~isempty(findstr(d2(d2i).name,'trypsindisp')) & ... 			% name has to contain 'trypsindisp'
       ~isempty(findstr(d2(d2i).name,'NoSine')), 			% name should contain 'NoSine'
      		d2b(si+2)	=d2(d2i);
		    d2b(si+2).name 	= d2(d2i).name;
            
           	%second(si) 	= str2num(name(18:length(name))) % many pairs, we should activate this line.
            si = si + 1;
      end;
   end;
    
   %====================================================================================

 % LOOP OVER FILES

for kloop = 3: length(d2b),
    eval(['cd C:\C2C12\Displacement\Syncronized\' cellfolder]); 
    displ = xlsread(d2b(kloop).name);
    savefilename = d2b(kloop).name; %from before, needed for saving figures
        
    m = findstr(savefilename,'.xls');
    savenumber = savefilename(m-5:m-1);
    time = str2num(savenumber)
%============================================================================================
% Find displacement values corresponding to position closest to magnetic bead
% Write displacement entries into a new excel file

I1 = find(displ(:,1) == x(1) & displ(:,2) == y(1));
data1(kloop-2,:) =  {num2str(time/10), x(1),y(1), ...
                        displ(I1,3), ...
                        displ(I1,4)};

I2 = find(displ(:,1) == x(2) & displ(:,2) == y(2));
data2(kloop-2,:) =  {num2str(time/10), x(2),y(2), ...
                        displ(I2,3), ...
                        displ(I2,4)};

I3 = find(displ(:,1) == x(3) & displ(:,2) == y(3));
data3(kloop-2,:) =  {num2str(time/10), x(3),y(3), ...
                        displ(I3,3), ...
                        displ(I3,4)};
                    
I4 = find(displ(:,1) == x(4) & displ(:,2) == y(4));
data4(kloop-2,:) =  {num2str(time/10), x(4),y(4), ...
                        displ(I4,3), ...
                        displ(I4,4)};
                    

end; %kloop



header   = '';       
 
     colnames = {'Time';'X';'Y'; 'DisplacementX'; ...
         'DisplacementY'};

     xlsfilename1 = 'Disparoundbead1.xls';
     xlsfilename2 = 'Disparoundbead2.xls';
     xlsfilename3 = 'Disparoundbead3.xls';
     xlsfilename4 = 'Disparoundmaxpt.xls';
     sheetname    = strcat(cellfolder,'_NoSine');
     
     cd('C:\C2C12\traction\Syncronized\');
     xlswrite(data1,header,colnames,xlsfilename1,sheetname);     
     xlswrite(data2,header,colnames,xlsfilename2,sheetname); 
     xlswrite(data3,header,colnames,xlsfilename3,sheetname); 
     xlswrite(data4,header,colnames,xlsfilename4,sheetname); 


end; %fil