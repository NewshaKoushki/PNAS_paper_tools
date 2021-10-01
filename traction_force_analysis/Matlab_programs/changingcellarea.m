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

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
fdate = 'croppeddata';
folds  = {'Cell33';'Cell34';'Cell76';'Cell67';'Cell78';'Cell44';'Cell63';'Cell77'};
%folds = {'Cell1','Cell2','Cell3','Cell4','Cell5','Cell6','Cell8'};
filename = {'cellboundary1';'cellboundary2';'cellboundary3'; 'cellboundary4'};
bdname = {'Magbeadpos1';'Magbeadpos2';'Magbeadpos3';'Magbeadpos4'};
ptname = {'selectpt1pos1';'Selectpt1pos2';'Selectpt1pos3';'Selectpt1pos4'};

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

for fi1 = 1 : length(folds),
   
  
   cellfolder	= folds{fi1}
%*   eval(['cd C:\C2C12\' fdate '\' cellfolder]); pwd,   
   eval(['cd C:\stretch\' fdate '\' cellfolder]); pwd,   
%*   eval(['Oct_18_01secondfiles']);         
   eval(['stretchinput']);         
   
      ButtonName1=questdlg(strcat('Do you want to draw cell boundary for ',cellfolder,' ?'), ...
                       'x,y coordinates of cell boundary','Yes');

   switch ButtonName1,
     case 'Yes',
         
   for bound = 1: length(cell_img),
        
       if bound == 1
       bound1=bound+1;
           else
            bound1=bound;
        end 
        
        if cell_img{bound} == 'crop_phase.tif',
            img = cell_img{bound};
            else
            img = strcat('cropCell1_',cell_img{bound});
        end

        
 % Use same cell boundary (drawn once) for cases with same cell boundary image
 
 if bound == 1 | ~strcmp(cell_img{bound1},cell_img{bound1-1})
   
     
   iptsetpref('truesize','off');
   figure(10); clf; hold on;
%   set(gcf,'Position',[1 29 1024 721],'NumberTitle','off','Name',pwd);
   set(gcf,'Position',[1 1 768 768],'NumberTitle','off','Name',pwd);
   set(gca,'Position',[0 0 1 1]);
   zoom on; pixval on; axis ij; 

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
   namefile=[filename{bound}]
   save(namefile,'xrub','yrub') 
 
     else
    % Copy the previous file to the new file
     copyfile(strcat(namefile,'.mat'),strcat(filename{bound},'.mat'))
   
%     close all;
%     clear all;
 

end; % if loop for bound    
end; % bound
end; % switch

% Select the x and y coordinates of the magnetic bead

   ButtonName=questdlg(strcat('Do you want to specify mag bead position for ',cellfolder,' ?'), ...
                       'x,y coordinates of magnetic bead','Yes');

   switch ButtonName,
     case 'Yes',
         
  for bound = 1: length(cell_img),
        
       if bound == 1
       bound1=bound+1;
           else
            bound1=bound;
        end 
        
        if cell_img{bound} == 'crop_phase.tif',
            img = cell_img{bound};
            else
            img = strcat('cropCell1_',cell_img{bound});
        end

 % Use same magnetic bead position(determined once) for cases with same cell boundary image
 
 if bound == 1 | ~strcmp(cell_img{bound1},cell_img{bound1-1})
   
         
   iptsetpref('truesize','off');
   figure(11); clf; hold on;
%   set(gcf,'Position',[1 29 1024 721],'NumberTitle','off','Name',pwd);
   set(gcf,'Position',[1 1 768 768],'NumberTitle','off','Name',pwd);
   set(gca,'Position',[0 0 1 1]);
   zoom on; pixval on; axis ij; 

   im = imread(img);
   im = double(im);
   im = im - min2(im);
   im = im ./ max(max(im)).*255;
   im = uint8(im);
   imshow(im);
   [BW, xrub1, yrub1] = roipoly(im);
   plot(xrub1,yrub1,'g+-'); 
   
   xt=num2str(xrub1);
   yt=num2str(yrub1);
   namefile=[bdname{bound}]
   save(namefile,'xrub1','yrub1') 

        else
    % Copy the previous file to the new file
     copyfile(strcat(namefile,'.mat'),strcat(bdname{bound},'.mat'))
   
     close all;
%     clear all;
 

end; % if loop for bound    
end; % bound

end; % switch

% Select the x and y coordinates of the point where tractions are to be determined

   ButtonName=questdlg(strcat('Do you want to specify another position for ',cellfolder,' ?'), ...
                       'x,y coordinates','Yes');

   switch ButtonName,
     case 'Yes',
         
  for bound = 1: length(cell_img),
        
       if bound == 1
       bound1=bound+1;
           else
            bound1=bound;
        end 
        
        if cell_img{bound} == 'crop_phase.tif',
            img = cell_img{bound};
            else
            img = strcat('cropCell1_',cell_img{bound});
        end

 % Use same point (determined once) for cases with same cell boundary image
 
 if bound == 1 | ~strcmp(cell_img{bound1},cell_img{bound1-1})
   
         
   iptsetpref('truesize','off');
   figure(12); clf; hold on;
   set(gcf,'Position',[1 1 768 768],'NumberTitle','off','Name',pwd);
%   set(gcf,'Position',[1 29 1024 721],'NumberTitle','off','Name',pwd);
   set(gca,'Position',[0 0 1 1]);
   zoom on; pixval on; axis ij; 

   im = imread(img);
   im = double(im);
   im = im - min2(im);
   im = im ./ max(max(im)).*255;
   im = uint8(im);
   imshow(im);
   [BW, xrub2, yrub2] = roipoly(im);
   plot(xrub2,yrub2,'g+-'); 
   
   xt=num2str(xrub2);
   yt=num2str(yrub2);
   namefile=[ptname{bound}]
   save(namefile,'xrub2','yrub2') 

        else
    % Copy the previous file to the new file
     copyfile(strcat(namefile,'.mat'),strcat(ptname{bound},'.mat'))
   
     close all;
%     clear all;
 

end; % if loop for bound    
end; % bound

end; % switch


end; % folds

