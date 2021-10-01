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
folds = {'Cell1002'};

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
    if  ~isempty(findstr(d2tryp(count).name,'Pos')),	% name must  contain Pos
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
   for d2i = 3 : length(d2),
      if ~isempty(findstr(d2(d2i).name,'.tif')) &...% name has to contain 'tif'
              ~isempty(findstr(d2(d2i).name,'t00')),
      		d2b(si+2)	=d2(d2i);
		                	si = si + 1;
        end;
   end;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Compute number of time points
timepts = (length(d2)-2)/(2*pos)

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    
% Draw boundary for cells as a function of position
     
   for bound = 3: length(d2b),
          
   figure(1); clf; hold on;
   set(gcf,'Position',[1 29 1024 721],'NumberTitle','off','Name',pwd);
    set(gcf,'Position',[1 1 768 768],'NumberTitle','off','Name',pwd);
    set(gca,'Position',[0 0 1 1]);
%    pixval on; axis ij; 
   phasetemp = d2b(bound).name;
   n=findstr(phasetemp,'_ch');
   phase = strcat('C:\SingleCell\croppeddata\',cellfolder,'\BigFB\',phasetemp(1:n-1),'_ch01.tif');

  % for cells = 1 : 100,
%           ButtonName1=questdlg(strcat('Do you want to draw cell boundaries for ',phasetemp(1:n-1),' ?'), ...
%                         'x,y coordinates of cell boundary','Yes');
%  
%     switch ButtonName1,
%       case 'Yes',
           img=phase;
           im = imread(img);
           im = double(im);
           im = im - min2(im);
           im = im ./ max(max(im)).*255;
           im = uint8(im);
           imshow(im);
%            [BW, xrub, yrub] = roipoly(im);
%            plot(xrub,yrub,'r.-'); 
% 
%            xt=num2str(xrub);
%            yt=num2str(yrub);
% 
%            namefile=strcat(phasetemp(1:n-1),'_',num2str(cells),'.mat')
%            save(namefile,'xrub','yrub') 

            [BW, xrub1, yrub1] = roipoly(im);
            plot(xrub1,yrub1,'b.-'); 
 
            xt=num2str(xrub1);
            yt=num2str(yrub1);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(2),'.mat')
            save(namefile,'xrub1','yrub1') ;
            
             [BW, xrub2, yrub2] = roipoly(im);
            plot(xrub2,yrub2,'g.-'); 
 
            xt=num2str(xrub2);
            yt=num2str(yrub2);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(3),'.mat')
            save(namefile,'xrub2','yrub2') 
            
             [BW, xrub3, yrub3] = roipoly(im);
            plot(xrub3,yrub3,'r.-'); 
 
            xt=num2str(xrub3);
            yt=num2str(yrub3);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(4),'.mat')
            save(namefile,'xrub3','yrub3')
            
             [BW, xrub4, yrub4] = roipoly(im);
            plot(xrub4,yrub4,'g.-'); 
 
            xt=num2str(xrub4);
            yt=num2str(yrub4);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(5),'.mat')
            save(namefile,'xrub4','yrub4')
            
             [BW, xrub5, yrub5] = roipoly(im);
            plot(xrub5,yrub5,'g.-'); 
 
            xt=num2str(xrub5);
            yt=num2str(yrub5);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(6),'.mat')
            save(namefile,'xrub5','yrub5')
            
             [BW, xrub6, yrub6] = roipoly(im);
            plot(xrub6,yrub6,'g.-'); 
 
            xt=num2str(xrub6);
            yt=num2str(yrub6);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(7),'.mat')
            save(namefile,'xrub6','yrub6')
            
             [BW, xrub7, yrub7] = roipoly(im);
            plot(xrub7,yrub7,'g.-'); 
 
            xt=num2str(xrub7);
            yt=num2str(yrub7);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(8),'.mat')
            save(namefile,'xrub7','yrub7')
            
             [BW, xrub8, yrub8] = roipoly(im);
            plot(xrub8,yrub8,'g.-'); 
 
            xt=num2str(xrub8);
            yt=num2str(yrub8);
 
            namefile=strcat(phasetemp(1:n-1),'_',num2str(9),'.mat')
            save(namefile,'xrub8','yrub8')
            
            

%            
%         case 'No'
%             break;
%     end;
   end;
   
end; 



% end; 

     close all;
%     clear all;

%end; % folds

