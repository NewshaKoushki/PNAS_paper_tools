%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the synchronous displacement profile
% Author: Ramaswamy Krishnan on 10/25/06
% Modified: 10/31/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

%====================================================================================
%  INPUT: LOAD DISPLACEMENTS

   close all;
   clear all;
   
     folds = {'Cell24'};
%     folds = {'Cell55'};
     matrix = [];
     
for fi1 = 1 : length(folds),
    cellfolder	= folds{fi1}
    eval(['Oct_18_01secondfiles']); 
    eval(['cd C:\C2C12\Displacement\Syncronized\' cellfolder]);
    pwd;
    savedirectory = strcat('C:\C2C12\traction\Syncronized\',cellfolder,'\');

%====================================================================================
% Settings for an animated movie

     mvfile3 = strcat(savedirectory,'Displacement',cellfolder,'.avi');
     aviobj3 = avifile(mvfile3,'compression','indeo5');

 aviobj3.Fps = 1;
 
 %============================================================================================
% Identify displacement files
   d2 = dir; si = 1;
   d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
   d2b(2)=d2(2);
   for d2i = 3 : length(d2),
      if ~isempty(findstr(d2(d2i).name,'.xls')) & ...		% search for an excel file
          ~isempty(findstr(d2(d2i).name,'trypsindisp')) & ... 			% name has to contain 'trypsindisp'
       isempty(findstr(d2(d2i).name,'NoSine')), 			% name should not contain 'NoSine'
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
    
%============================================================================================



%  CROP DISPLACEMENTS INTO SIZE (even N)^2

sizx  	= 	size(unique(displ(:,1)),1);
sizy  	= 	size(unique(displ(:,2)),1);
minsiz 	= 	min(sizx,sizy);
if mod(minsiz,2) ~= 0, minsiz = minsiz - 1; end;
displ2 	= 	zeros(minsiz^2,4);
for i = 1 : 4,
	displ1 = reshape(displ(:,i),sizy,sizx);
	displ1 = displ1(1 : minsiz,1 : minsiz);
	displ2(:,i) = reshape(displ1,minsiz^2,1);
end;
displ 	= 	displ2;

%====================================================================================
%  CENTER DISPLACEMENTS, PIXELS TO MICRONS, SPACING

xv  = (displ(:,1) - mean(displ(:,1))) * pixelsize;  
yv  = (displ(:,2) - mean(displ(:,2))) * pixelsize;
uxv = (displ(:,3)) * pixelsize;
uyv = (displ(:,4)) * pixelsize;

xvorig = xv;
yvorig = yv;

spacing = yv(2) - yv(1);

%====================================================================================
%  N, RESHAPE DISPLACEMENTS

N   = sqrt(length(xv));

 x  = reshape( xv,N,N);
 y  = reshape( yv,N,N);
ux  = reshape(uxv,N,N);
uy  = reshape(uyv,N,N);


%=====================================================================================
%Filter for the displacement matrix

ux  = myfilter_exp(ux, 15); % f0 = 15, the cut-off frequency.
uy  = myfilter_exp(uy, 15);

%====================================================================================
%  Save figures parameters  


     m = findstr(savefilename,'.xls');
     savenumber = savefilename(m-5:m-1);
     time = str2num(savenumber)
     
     for freqloop = 1: length(freq)
         
         if time >= min_time(freqloop) & time <= max_time(freqloop)
             frequency = freq(freqloop);
         end %(if loop)
     end %(freqloop)
   
     disptitle = strcat('Displacement;',' time = ',num2str(time/10),'s; freq = ',num2str(frequency));
%    
%      disptitle='Displacement';

%====================================================================================
%  CELL BOUNDARY
%      start= (findstr(d2b(k).name,'Cell'));		% name has to contains 'Cell'
%      ed=    findstr(d2b(k).name,'disp');		% name has to cotains 'disp'
%             m=d2b(k).name(start+4:ed-1);
% where m is the number of the cell folder

 eval(['cd C:\C2C12\croppeddata\' cellfolder ]); pwd;
    
    
    
 for cellbound = 1: length(min_time)
 
     if str2num(savenumber) >= min_time(cellbound) & str2num(savenumber) <= max_time(cellbound)
     filename=strcat('cellboundary',num2str(cellbound),'.mat');
     filename1 = strcat('Magbeadpos',num2str(cellbound),'.mat');
     filename2 = strcat('Selectpt1pos',num2str(cellbound),'.mat');
 %%[filename,pathname] = uigetfile('*.*','Cell boundary');
 %%  cd(pathname);
   load(filename);
   load(filename1);
   load(filename2);
   xrub = (xrub - mean(displ(:,1))) * pixelsize; 
   yrub = (yrub - mean(displ(:,2))) * pixelsize;
   xrub1 = (xrub1 - mean(displ(:,1))) * pixelsize; 
   yrub1 = (yrub1 - mean(displ(:,2))) * pixelsize;
   xrub2 = (xrub2 - mean(displ(:,1))) * pixelsize; 
   yrub2 = (yrub2 - mean(displ(:,2))) * pixelsize;

   frequency = freq(cellbound);
end %(if loop)
end %(cellbound)
   
area_cell = polyarea(xrub,yrub);

cell1 = find(xv >= min(xvorig) & xv <= max(xvorig) &...
             yv >= min(yvorig) & yv <= max(yvorig));
cell  = inpolygon(xv(cell1),yv(cell1),xrub,yrub);
cell  = cell1(find(cell > 0));       % POINTS IN THE INTERIOR OF THE CELL

%====================================================================================

% Save displacement figures and movie
       figure(1); 
   rect = [100, 100, 800, 560]; 
   set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;
%    CAXIS([0 0.12])

   surf(x,y,sqrt(ux.^2+uy.^2)); 
   view(2); colormap jet; shading interp; 
   cbh = colorbar; set(cbh,'FontSize',12,'FontWeight','bold');
   m1 = max2(sqrt(ux.^2+uy.^2));
   h = quiver3(x,y,m1*ones(size(x)),ux,uy,zeros(size(x)),1);
   set(h,'Color','w','LineWidth',1);
   axis equal; axis tight; axis ij; 
%  axis([-30 30 -30 30]);    
   axis([min2(x) max2(x) min2(y) max2(y)])
   xlabel('x (\mum)','FontSize',12,'FontWeight','bold'); ylabel('y (\mum)','FontSize',12,'FontWeight','bold');
   set(gca,'FontSize',12,'FontWeight','bold'); box on;
   title(disptitle,'FontSize',12,'FontWeight','bold'); 
   
   % Save displacement figures and movie
       figure(2); 
   rect = [100, 100, 800, 560]; 
   set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;
%   CAXIS([0 0.5])
   surf(x,y,sqrt(ux.^2+uy.^2)); 
   view(2); colormap jet; shading interp; 
   cbh = colorbar; set(cbh,'FontSize',12,'FontWeight','bold');
   m1 = max2(sqrt(ux.^2+uy.^2));
   h = quiver3(x,y,m1*ones(size(x)),ux,uy,zeros(size(x)),1);
   set(h,'Color','w','LineWidth',1);
   axis equal; axis tight; axis ij; 
%  axis([-30 30 -30 30]);    
   axis([min2(x) max2(x) min2(y) max2(y)])
   xlabel('x (\mum)','FontSize',12,'FontWeight','bold'); ylabel('y (\mum)','FontSize',12,'FontWeight','bold');
   set(gca,'FontSize',12,'FontWeight','bold'); box on;
   title(disptitle,'FontSize',12,'FontWeight','bold'); 
   
%====================================================================================
%  Save figures parameters  
   
%    m = findstr(savefilename,'disp');
%    n = findstr(savefilename,'.');
%    savedate = savefilename(1:m-1);
%    savedate(findstr(savedate,'-')) = '_';
%    savenumber = savefilename(m+4:n-1);
% %   disptitle = strcat('Displacement',savenumber);
% %   tractitle = strcat('Constrained FTTC: Tractions',savenumber);
%    
%    
% %====================================================================================
% %  SAVE FIGURES
%    m = findstr(savefilename,'disp');
%    n = findstr(savefilename,'.');
%    savedate = savefilename(1:m-1);
%    savedate(findstr(savedate,'-')) = '_';
%    savenumber = savefilename(m+4:n-1);
%     saveas(1,strcat(savedirectory,savedate,'_Disp',savenumber,'.fig'));
% %    saveas(2,strcat(savedirectory,savedate,'_Unconstraint_Traction',savenumber,'.fig'));   
% %    saveas(3,strcat(savedirectory,savedate,'_Constraint_Traction',savenumber,'.fig'));
% %    saveas(4,strcat(savedirectory,savedate,'_Result',savenumber,'.fig'));

%====================================================================================
   
    figure(1);
   h = plot3(xrub,yrub,m1*ones(size(xrub)),'w-'); set(h,'LineWidth',2.0);

   aviobj3 = addframe(aviobj3,figure(1));  
   
end; %(kloop)
   
   aviobj3 = close(aviobj3);
   
   clear a* A* B* c* C* d* D* h* i* j* k* m* n* p* r* s* T* t* u* U* v* x* y*;
 close all;

end; %fil

return;