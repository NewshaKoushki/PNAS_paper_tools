%DRAW_BOUNDARY is used for selecting the vertices of a polygon
%	that constitutes the cell boundary.
%	Use normal button clicks to add vertices to the polygon. 
%	Pressing <BACKSPACE> or <DELETE> removes the previously selected
%	vertex. A shift-click, right-click, or double-click adds a
%	final vertex to the selection and then starts the fill;
%	pressing <RETURN> finishes the selection without adding a
%	vertex.
%	The x- and y-coordinates (xrub and yrub) of the vertices are 
%	saved in a file called cellboundary in the current folder.

%	Iva Marija Tolic-Norrelykke 03-21-01

clear all;
close all;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
fdate = 'croppeddata';
% folds = {'Cell145';'Cell146';'Cell147';'Cell150';'Cell151';'Cell152'};
folds = {'Cell0518'};

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

for fi1 = 1 : length(folds),
    
  
   cellfolder	= folds{fi1}
   
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   % No. of positions
nchannels = 3;

eval(['cd C:\SingleCell\rawdata' '\' cellfolder, '\Mark_and_Find_002']); 
    d2tryp = dir;
    si = 0,
    for count = 1 : length(d2tryp),
    if  ~isempty(findstr(d2tryp(count).name,'Position')),	% name must  contain Pos
          	si = si + 1;
      end;    
    end;
    pos = si/nchannels
 %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   eval(['cd C:\SingleCell\displacement\' cellfolder ]); pwd,
   
%    eval(['cd C:\SingleCell\croppeddata\' cellfolder '\BF\']); pwd,   
%    %eval(['stretchinput']);         
%    
%       ButtonName1=questdlg(strcat('Do you want to draw cell boundaries for ',cellfolder,' ?'), ...
%                        'x,y coordinates of cell boundary','Yes');
% 
%    switch ButtonName1,
%      case 'Yes',
         
   %  eval(['cd C:\uniaxialstretch\croppeddata'  '\' cellfolder  ]); pwd;
     d2 = dir; si = 1;
   d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
   d2b(2)=d2(2);
  % this will be used to define number of time points
   for d2i = 3 : length(d2),
      if ~isempty(findstr(d2(d2i).name,'.tif'))% &...% name has to contain 'tif'
              ~isempty(findstr(d2(d2i).name,'t00')),
      		d2b(si+2)	=d2(d2i);
		                	si = si + 1;
        end;
   end;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Compute number of time points
timepts = (length(d2)-2)/(2*pos)
%timepts=(si-2);  % this is what I (Hossein) thinks should be here.
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    
% Draw boundary for cells as a function of position
     
 for bound = 3: length(d2b),
     
   if bound==3     
   figure(1); clf; hold on;
   set(gcf,'Position',[1 29 1024 721],'NumberTitle','off','Name',pwd);
    set(gcf,'Position',[1 1 768 768],'NumberTitle','off','Name',pwd);
    set(gca,'Position',[0 0 1 1]);
%    pixval on; axis ij; 
   phasetemp = d2b(bound).name;
   n=findstr(phasetemp,'_ch');
   phase = strcat('C:\SingleCell\croppeddata\',cellfolder,'\BigFB\',phasetemp(1:n-1),'_ch02.tif');
%%%%%%%%
    thenumberofcells=0% hossein added this line
       for cells = 1 : 100,
          ButtonName1=questdlg(strcat('Do you want to draw cell boundaries for ',phasetemp(1:n-1),' ?'), ...
                        'x,y coordinates of cell boundary','Yes');
     switch ButtonName1,
      case 'Yes',
           img=phase;
           im = imread(img);
           im = double(im);
           im = im - min2(im);
           im = im ./ max(max(im)).*255;
           im = uint8(im);
           imshow(im);
           [BW, xrub, yrub] = roipoly(im);
           plot(xrub,yrub,'r.-'); 

           xt=num2str(xrub);
           yt=num2str(yrub);

           namefile=strcat(phasetemp(1:n-1),'_',num2str(cells),'.mat')
           save(namefile,'xrub','yrub') 
           thenumberofcells=thenumberofcells+1
        case 'No'
            break;
    end;
   end;
   else  % hossein added this line
           
        for i=1:thenumberofcells
            A1=load('05182015_TF_9kPa_Position001_t0000_1.mat')
            A2=load('05182015_TF_9kPa_Position001_t0000_2.mat')
            A3=load('05182015_TF_9kPa_Position001_t0000_3.mat')
            A4=load('05182015_TF_9kPa_Position001_t0000_4.mat')
            A5=load('05182015_TF_9kPa_Position001_t0000_5.mat')
            A6=load('05182015_TF_9kPa_Position001_t0000_6.mat')
               time=bound-3
                filename1=horzcat('05182015_TF_9kPa_Position001_t',num2str(time,'%04d'),'_1.mat');
                filename2=horzcat('05182015_TF_9kPa_Position001_t',num2str(time,'%04d'),'_2.mat');
                filename3=horzcat('05182015_TF_9kPa_Position001_t',num2str(time,'%04d'),'_3.mat');
                filename4=horzcat('05182015_TF_9kPa_Position001_t',num2str(time,'%04d'),'_4.mat');
                filename5=horzcat('05182015_TF_9kPa_Position001_t',num2str(time,'%04d'),'_5.mat');
                filename6=horzcat('05182015_TF_9kPa_Position001_t',num2str(time,'%04d'),'_6.mat');
                save(filename1,'A1');
                save(filename2,'A2');
                save(filename3,'A3');
                save(filename4,'A4');
                save(filename5,'A5');
                save(filename6,'A6');
            end
        
        
        
        
        
        
   end
   
   
   
   
   
   
   
end; 

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Reproduce cell boundary at time t=0 for all time points for a given position
  eval(['cd C:\SingleCell\displacement\' cellfolder ]); pwd, 
     d2bound = dir; sibound = 1;
   d2bbound(1)=d2bound(1);			%d2b is part of d2 and will include [.][..] and original cell files
   d2bbound(2)=d2bound(2);
   for d2ibound = 3 : length(d2bound),
      if ~isempty(findstr(d2bound(d2ibound).name,'.mat')) &...% name has to contain 'mat'
              ~isempty(findstr(d2bound(d2ibound).name,'t00')),
      		d2bbound(sibound+2)	=d2bound(d2ibound);
		                	sibound = sibound + 1;
        end;
   end;
%%%%%%% Hossein: what is this section doing???????
   for i = 3 : length(d2bbound),
       for j = 1 : timepts-1,
        if j <=9,
           copyfile(d2bbound(i).name, strcat(d2bbound(i).name(1:31),num2str(j), d2bbound(i).name(33:38)));
        else
            copyfile(d2bbound(i).name, strcat(d2bbound(i).name(1:30),num2str(j), d2bbound(i).name(33:38)));
        end;
       end;
   end;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    

end; 

     close all;
%     clear all;

%end; % folds

