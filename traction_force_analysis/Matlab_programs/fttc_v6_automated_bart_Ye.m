%  FTTC calculates the tractions, using both the Unconstrained
%	and the Constrained FTTC.
%	The required input is a file of displacements (4 columns:
%	x-coordinates,
%	y-coordinates, x-displacements, y-displacements, all in pixels), 
%	pixel to micron conversion factor, Young's modulus and Poisson's ratio 
%	of the gel, and the cell boundary file (it is needed for the 
%	Constrained FTTC only).

%	Iva Marija Tolic-Norrelykke 03-22-01; edited by Ramaswamy Krishnan
%	06-26-12
   
% Trac o/p: Filename; RMS traction; Median traction; CM; Cell area 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


   close all;
   clear all;
     
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%  INPUT: PIXELSIZE, YOUNG'S MODULUS, POISSON'S RATIO


folds         = {'Cell157'};
% pixelsize     = 0.462; %0.92 for Rama's Inverted 10x ; 1.297 for Xinyong 10x; 0.462 for Bart - 20x; 0.2522 for Allen dynamic; 0.3027 for Allen Static
%  young        = 26000; %Pa
%  pois         = 0.48;
%  ncells       = 1;
%  young         = young * 1e-12; % do not change
% RMSTmin       = 10; % don't change

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
 

    cellfolder	= folds{1}
    
      eval(['cropinput']);
      
        pois         = 0.48;
  ncells       = 1;
  young         = young * 1e-12; % do not change
 RMSTmin       = 10; % don't change
      
eval(['cd C:\Singlecell\traction\' cellfolder]); 
      str1 = 'Pos     Side    t       Cell    RMST    MT      CM            F    Area  L(1,1) L(2,2) Theta';
     fName = 'Results.txt';         %# A file name
fid = fopen(fName,'w');            %# Open the file
 if fid ~= -1
 %  fprintf(fid,'%s\r\n',str);       %# Print the string
   fprintf(fid,'%s\r\n',str1);       %# Print the string
   fclose(fid);                     %# Close the file
 end

 eval(['cd C:\Singlecell\displacement\' cellfolder]); 
     savedirectory = strcat('C:\Singlecell\traction\',cellfolder,'\');
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Settings for an animated movie

     mvfile1 = strcat(savedirectory,'Constraint_Traction',cellfolder,'_',num2str(1),'.avi');
     mvfile2 = strcat(savedirectory,'UnConstr_Traction',cellfolder,'_',num2str(1),'.avi');
     mvfile3 = strcat(savedirectory,'Displacement',cellfolder,'_',num2str(1),'.avi');
     aviobj1 = avifile(mvfile1,'compression','None');
     aviobj2 = avifile(mvfile2,'compression','None');
     aviobj3 = avifile(mvfile3,'compression','None');

      aviobj1.Fps = 2;
      aviobj2.Fps = 2;
      aviobj3.Fps = 2;
%   
   
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
 
 

%====================================================================================
%  CONSTANTS

a1   = (1.0 + pois) * (1.0 - pois) / (pi * young);
b1   = (1.0 + pois) * pois / (pi * young);
c1   = (1.0 + pois) * pois / (pi * young);

% LOOP OVER FILES

for kloop = 3: length(d2b),
    eval(['cd C:\Singlecell\displacement\' cellfolder]); 
    displ = load(d2b(kloop).name);
    savefilename = d2b(kloop).name %from before, needed for saving figures
   

%====================================================================================

% To compute number of cells per displacement image (corresponds to number
% % of cell boundary files in the ../Displacements folder for each image
%    d2disp = dir; si = 1; ncells = 0;
%    d2bdisp(1)=d2disp(1);			%d2b is part of d2 and will include [.][..] and original cell files
%    d2bdisp(2)=d2disp(2);
%    for d2idisp = 3 : length(d2disp),
%       if ~isempty(findstr(d2disp(d2idisp).name,'.mat')) &  ...		% search for a mat file
%           ~isempty(findstr(d2disp(d2idisp).name,savefilename(1:15))), 			% name has to contain particular string
%       		d2bdisp(si+2)	=d2disp(d2i);
% 		    d2bdisp(si+2).name 	= d2disp(d2i).name;
%             si = si + 1;
%             ncells = ncells + 1;
%       end;
%    end;
%  
  digit = findstr(savefilename, '_ch');
    
  
%====================================================================================

        for   cellno = 1 : ncells,
            
           
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
%  xv  = (displ(:,1)) * pixelsize;  
%  yv  = (displ(:,2)) * pixelsize;

 
uxv = (displ(:,3)) * pixelsize;
uyv = (displ(:,4)) * pixelsize;

xvorig = xv;
yvorig = yv;

spacing = yv(2) - yv(1);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%     UNCONSTRAINED FTTC
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
%  kx, ky, AND  k_abs 

clear mala*;
for i = 1 : (N/2)
   malax(i,:) = 0 : ((N/2)-1);                     
   malay(i,1:(N/2)) = (N/2)-i;                   
end;

kx = [ malax  malax-(N/2);  malax  malax-(N/2) ];
ky = [ malay-(N/2)  malay-(N/2);  malay  malay ];
ky = flipud(ky);
kx(:,(N/2+1)) =  kx(:,(N/2+1));                     
ky((N/2+1),:) =  ky((N/2+1),:);

k_abs = sqrt(kx.^2 + ky.^2);

%====================================================================================
%  ALPHA

alpha = atan2(ky,kx);
if kx(1,1) == 0 & ky(1,1) == 0,
   alpha(1,1) = 1.57080;
end;

%====================================================================================
%  C AND D

Cx = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (cos(alpha)).^2);
Cy = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (sin(alpha)).^2);
D  = ((k_abs * young) / (2 * (1 - pois^2))) .* (pois .* sin(alpha) .* cos(alpha));

D(:,(N/2+1)) = zeros(N,1);
D((N/2+1),:) = zeros(1,N);

%====================================================================================
%  CALCULATE THE TRACTIONS

Dx = fft2(ux * 2 * pi / (N * spacing));
Dy = fft2(uy * 2 * pi / (N * spacing));

Tx = Cx .* Dx  +  D  .* Dy;
Ty = D  .* Dx  +  Cy .* Dy;

tx = real(ifft2(Tx));
ty = real(ifft2(Ty));
     
%====================================================================================
% RMS TRACTION, THETA, TRACE, U

rmst_onestep = sqrt(mean2(tx.^2 + ty.^2));

% average over the large region = sum_over_region(fraction_i* force_i)/sum_over_region(force_i)

clear i;
Dxx = imag(Tx(1,2)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6; 
Dyy = imag(Ty(2,1)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6; 
Dxy = imag(Ty(1,2)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6; 
Dyx = imag(Tx(2,1)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6; 

[a_sym, lambda, theta] = rotate1([Dxx Dxy; Dyx Dyy]);

Trace_moment = Dxx + Dyy;

Uonestep = 0.5 * sum(sum([tx ty] .* [ux uy])) * spacing^2 * 1e-6;

interior = find(x ~= min2(x) & x ~= max2(x) & y ~= min2(y) & y ~= max2(y));
Uonestepint = 0.5 * sum(sum([tx(interior) ty(interior)].*[ux(interior) uy(interior)])) * spacing^2 * 1e-6;

%====================================================================================
%  Save figures parameters  


     m = findstr(savefilename,'.xls');
     disptitle = strcat('Displacement(Microns)');
     tractitle = strcat('Traction map (Pa)');
     untractitle = strcat('Unconstrained FTTC: Tractions (Pa)');
   
%====================================================================================
%  FIGURES

   figure(1); 
   
     rect = [100, 100, 800, 560]; 
     set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;

  % caxis([0 0.3])
   surf(x,y,sqrt(ux.^2+uy.^2)); 
   view(2); shading interp; 

   cbh = colorbar; set(cbh,'FontSize',18,'FontWeight','bold');
   m1 = max2(sqrt(ux.^2+uy.^2))
   h = quiver3(x,y,m1*ones(size(x)),ux,uy,zeros(size(x)),1);
   set(h,'Color','w','LineWidth',1);
   axis equal; axis tight; axis ij; 
%  axis([-30 30 -30 30]);    
   axis([min2(x) max2(x) min2(y) max2(y)])
%   xlabel('x (\mum)','FontSize',18,'FontWeight','bold'); ylabel('y (\mum)','FontSize',18,'FontWeight','bold');
   set(gca,'FontSize',18,'FontWeight','bold'); box on;
    disptime = str2num(savefilename(14:15))*5;
    headerdisp = strcat(disptitle,';',savefilename(1:6),';',...
        num2str(disptime),'min');
   title(headerdisp,'FontSize',18,'FontWeight','bold');   
   
   
   figure(2); set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;
   %caxis([0 700])
   surf(x,y,sqrt(tx.^2+ty.^2)*1e12); 
   view(2); colormap jet; shading interp; 
   cbh = colorbar; set(cbh,'FontSize',18,'FontWeight','bold');
   m2 = max2(sqrt(tx.^2+ty.^2)*1e12);
   h = quiver3(x,y,m2*ones(size(x)),tx,ty,zeros(size(x)),1);
   set(h,'Color','w','LineWidth',1);
   axis equal; axis tight; axis ij;
   axis([min2(x) max2(x) min2(y) max2(y)]);
%   xlabel('x (\mum)','FontSize',18,'FontWeight','bold'); ylabel('y (\mum)','FontSize',18,'FontWeight','bold');
   set(gca,'FontSize',18,'FontWeight','bold'); box on;
   title(untractitle,'FontSize',18,'FontWeight','bold');
   
%    figure(4); set(gcf,'Position',rect);
%    clf; hold on; set(gca,'Visible','off');
%    left = 0.2; left2 = 0.5; mov = 0.05; top = 1; difv = 0.05;
%    text(left,top,'Pixel to \mum:','FontSize',18,'FontWeight','bold'); 
%    text(left2,top,num2str(pixelsize,'%10.3f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-difv,'Young''s modulus (Pa):','FontSize',18,'FontWeight','bold'); 
%    text(left2,top-difv,num2str(young*1e12,'%10.0f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-2*difv,'Poisson''s ratio:','FontSize',18,'FontWeight','bold'); 
%    text(left2,top-2*difv,num2str(pois,'%10.3f'),'FontSize',18,'FontWeight','bold');
%    left = -0.1; top = top-5*difv;
%    h = text(left+mov,top+difv,'Unconstrained FFTC'); set(h,'FontSize',18,'FontWeight','bold');
%    text(left,top-difv,'RMS traction (Pa):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-2*difv,num2str(rmst_onestep*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-3*difv,'Orientation of principle tractions (^o):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-4*difv,num2str(theta*180/pi,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-5*difv,'Net contractile moment (pJ):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-6*difv,num2str(-Trace_moment*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-7*difv,'Total strain energy (pJ):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-8*difv,num2str(Uonestepint*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
   
    
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%     CONSTRAINED FTTC
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%====================================================================================
%  N, RESHAPE DISPLACEMENTS

N = sqrt(length(xv)) * 2;

iteration = 100; 
the_same  = 1e-6;
   
ux    = ux * 2 * pi / (N * spacing);
uy    = uy * 2 * pi / (N * spacing);

[x,y]  = meshgrid([(min(xv) - spacing*N/4 : spacing : min(xv)), x(1,2:size(x,2)-1),...
         (max(xv) : spacing : max(xv) + spacing*N/4)],...
			[(min(yv) - spacing*N/4 : spacing : min(yv)), y(2:size(y,1)-1,1)',...
   		(max(yv) : spacing : max(yv) + spacing*N/4)]);
ux = [zeros(N/4,N); zeros(N/2,N/4), ux, zeros(N/2,N/4); zeros(N/4,N)];   
uy = [zeros(N/4,N); zeros(N/2,N/4), uy, zeros(N/2,N/4); zeros(N/4,N)];   
   
xv     =  reshape(x,N^2,1);
yv     =  reshape(y,N^2,1);
uxv    =  reshape(ux,N^2,1);
uyv    =  reshape(uy,N^2,1);

%====================================================================================
%  CELL BOUNDARY

eval(['cd C:\Singlecell\displacement\' cellfolder]); pwd;
      
   filename = strcat(savefilename(1:digit-1),'_',num2str(cellno),'.mat');
   load(filename);
   xrub = (xrub - mean(displ(:,1))) * pixelsize;
   yrub = (yrub - mean(displ(:,2))) * pixelsize;
 
area_cell = polyarea(xrub,yrub);

cell1 = find(xv >= min(xvorig) & xv <= max(xvorig) &...
             yv >= min(yvorig) & yv <= max(yvorig));
cell  = inpolygon(xv(cell1),yv(cell1),xrub,yrub);
cell  = cell1(find(cell > 0));       % POINTS IN THE INTERIOR OF THE CELL

%====================================================================================
%  kx, ky, AND k_abs 

clear mala*;
for i = 1:(N/2)
   malax(i,:) = 0:((N/2)-1);                     
   malay(i,1:(N/2)) = (N/2)-i;                   
end;

kx = [ malax  malax-(N/2);  malax  malax-(N/2) ];
ky = [ malay-(N/2)  malay-(N/2);  malay  malay ];
ky = flipud(ky);
kx(:,(N/2+1)) =  kx(:,(N/2+1));                     
ky((N/2+1),:) =  ky((N/2+1),:);

k_abs = sqrt(kx.^2 + ky.^2);

%====================================================================================
%  ALPHA

alpha = atan2(ky,kx);
if kx(1,1) == 0 & ky(1,1) == 0,
   alpha(1,1) = 1.57080;
end;

%====================================================================================
%  C AND D

Cx = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (cos(alpha)).^2);
Cy = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (sin(alpha)).^2);
D  = ((k_abs * young) / (2 * (1 - pois^2))) .* (pois .* sin(alpha) .* cos(alpha));

D(:,(N/2+1)) = zeros(N,1);
D((N/2+1),:) = zeros(1,N);

%====================================================================================
%  A AND B

if k_abs(1,1) == 0,
   k_abs(1,1) = 1;
end;

Ax = a1 * ( 2 * pi ./ k_abs) + b1 * (2 * pi ./ k_abs) .* (sin(alpha)).^2;
Ay = a1 * ( 2 * pi ./ k_abs) + b1 * (2 * pi ./ k_abs) .* (cos(alpha)).^2;
B  = c1 * (-2 * pi ./ k_abs) .* sin(alpha) .* cos(alpha);

B(:,(N/2+1)) = zeros(N,1);
B((N/2+1),:) = zeros(1,N);

%====================================================================================
%  CALCULATE THE FIRST TRACTIONS

Dx = fft2(ux);
Dy = fft2(uy);

Tx = Cx .* Dx  +  D  .* Dy;
Ty = D  .* Dx  +  Cy .* Dy;

tx = real(ifft2(Tx));
ty = real(ifft2(Ty));

%====================================================================================
%  ITERATIONS

step = 1;
are_the_same = 0;

while step <= iteration & are_the_same == 0,
   
%====================================================================================
%  FORM NEW (MIXED) TRACTION MATRIX

   tx2 = zeros(N,N);
   ty2 = zeros(N,N);

   tx2(cell) = tx(cell);   
   ty2(cell) = ty(cell);   
   
%====================================================================================
%  CALCULATE THE INDUCED DISPLACEMENTS

Tx = fft2(tx2);       
Ty = fft2(ty2);      

Dx = Ax .* Tx  +  B  .* Ty;  
Dy = B  .* Tx  +  Ay .* Ty;  

dx2 = real(ifft2(Dx));
dy2 = real(ifft2(Dy));

%====================================================================================
%  FORM NEW (MIXED) DISPLACEMENT MATRIX

   dx3 = dx2;                   
   dy3 = dy2;   

   dx3(cell) = ux(cell);  
   dy3(cell) = uy(cell);

%====================================================================================
%  CALCULATE THE TRACTIONS

Dx = fft2(dx3);
Dy = fft2(dy3);

Tx = Cx .* Dx  +  D  .* Dy;
Ty = D  .* Dx  +  Cy .* Dy;

tx3 = real(ifft2(Tx));
ty3 = real(ifft2(Ty));

%====================================================================================
%  FORM NEW (MIXED) TRACTION MATRIX
   
   tx = 0.5 * tx3 + 0.5 * tx2;
   ty = 0.5 * ty3 + 0.5 * ty2;
   
if spacing <= 4,
   tx = 0.2 * tx3 + 0.8 * tx2;
   ty = 0.2 * ty3 + 0.8 * ty2;
end;

maxt(step) = max(max(tx3.^2 + ty3.^2));

%====================================================================================
% CHECK IF THE TRACTIONS HAVE CONVERGED

if step > 3,
   if abs(maxt(step) - maxt(step-1)) <= the_same * maxt(step),
      are_the_same = 1; 
   end;
end;
step = step + 1;
end; %(iteration)
if are_the_same ~= 1, disp(' Convergence was not reached!'); end;

%====================================================================================
%  SUBTRACT mean(t(cell)), BACK TO ORIGINAL ux 
  
tx = real(tx3); 
ty = real(ty3);
tx(cell) = tx(cell) - mean(tx(cell));
ty(cell) = ty(cell) - mean(ty(cell));

   
tx4      = tx(cell);
ty4      = ty(cell);

   
   
ux = ux * N * spacing / (2 * pi);
uy = uy * N * spacing / (2 * pi);

% % Save tx(cell) and ty(cell) into a file
%         tracfilename = strcat('C:\C2C12\traction\',cellfolder,'\tracmap_',num2str(time));
%         fid = fopen(tracfilename,'wt');     % 'wt' means "write text"
%         if (fid < 0)
%             error('could not open file "tracmap_time.txt"');
%         end;
%         tracval = [sqrt(tx(cell).^2+ty(cell).^2)*1e12];
%         fprintf(fid,'%10.5f\n',tracval);
%         fclose(fid);
        
        
%====================================================================================
% %  RMS TRACTION, TOTAL FORCE, THETA, A_CUTS, TRACE, U


  count  = 1;
       for loop = 1: size(cell,1),   
          if or(abs(tx4(loop)*1e12) >=RMSTmin,abs(ty4(loop)*1e12) >=RMSTmin),
           tx5(count)        = tx4(loop);
           ty5(count)        = ty4(loop);
           count             = count + 1;
           data6(loop,:) =  {tx4(loop), ty4(loop)};
           
       end;
    end;
 
   
median_trac = median(sqrt(tx(cell).^2 + ty(cell).^2));   
rmst_iterative = sqrt(mean2(tx(cell).^2 + ty(cell).^2));      

% clear tx4,ty4,tx5,ty5


tot_force   = sum(sum(sqrt(tx(cell).^2 + ty(cell).^2)))*spacing^2;

Dxx = sum(sum([tx] .* [x])) * spacing^2 * 1e-6; 
Dxy = sum(sum([ty] .* [x])) * spacing^2 * 1e-6;
Dyx = sum(sum([tx] .* [y])) * spacing^2 * 1e-6;
Dyy = sum(sum([ty] .* [y])) * spacing^2 * 1e-6;

[a_sym, lambda, theta0] = rotate1([Dxx Dxy; Dyx Dyy]);

if abs(lambda(1,1)) > abs(lambda(2,2)), theta = -theta0;
else theta = pi/4 - theta0;
end;

a_tractions = tan(theta);
a_cuts = -1 / tan(theta);
   
Trace_moment = Dxx + Dyy;

Uecm = 0.5 * sum(sum([tx(cell) ty(cell)] .* [ux(cell) uy(cell)])) * spacing^2 * 1e-6;


%======================================================================================
%  CALCULATION OF PRESTRESS

max_height = 6;      % height of the cell at the nucleus
min_height = 0.5;    % height of the cell at the edges of the cell

%======================================================================================
%  POINTS OF CUTS ON THE LINE THROUGH (0,0)

number = 0;
clear point_x point_y

for i = floor(-sqrt(2)*N/2) : ceil(sqrt(2)*N/2),
   
   p = [i * spacing * cos(theta), i * spacing * sin(theta)];
   axbo = max([max(xrub),max(yrub),-min(xrub),-min(yrub)]);
   if inpolygon(p(1),p(2),[-axbo axbo axbo -axbo],...
         [-axbo -axbo axbo axbo]) > 0, 
   	number = number + 1;   
   	point_x(number) = p(1);
      point_y(number) = p(2);
   end; %(if)
   
end; %(for)

%======================================================================================
%  HEIGHTS IN POINTS

height = zeros(number,1);
middle = find(point_x == 0 & point_y == 0);
height(middle) = max_height;
for i = 1 : middle - 1,
   height(i) = min_height + (i-1) * (max_height - min_height) / (middle - 1);
end;
for i = number : (-1) : middle + 1,
   height(i) = min_height + (number-i) * (max_height - min_height) / (number - middle);
end;

%======================================================================================
%  PRESTRESS

clear sumforce diameter area_cross prestress

for i = 1 : number,
   
%===================================================================================
   %  B_CUTS, SUMFORCE
   
   b_cuts      = point_y(i) - a_cuts * point_x(i);
   if a_cuts > 0, znak = '>'; else znak = '<'; end; 
   eval(['halfplane = find(y ' znak ' a_cuts * x + b_cuts);']);
   cutleft     = intersect(halfplane,cell);
   if ~isempty(cutleft) & ~isempty(setdiff(cell,cutleft)),
   vec         = [sqrt(1/(1+a_tractions^2)) sqrt(1/(1+a_tractions^2))*a_tractions];
   clear vec_line; 
   for j = 1 : length(cutleft), 
      vec_line(j,:) = vec; 
   end;

   sumforce(i) = sum(sum([tx(cutleft)  ty(cutleft)]  .* vec_line)) * spacing^2 * 1e12;
   	%===================================================================================
   %  DIAMETER
   
   if abs(a_cuts) < 1,
   	xdots    = min(xrub) : 0.1 : max(xrub);
      ydots    = a_cuts .* xdots + b_cuts;
   else
   	ydots    = min(yrub) : 0.1 : max(yrub);
      xdots    = (ydots - b_cuts) ./ a_cuts;
   end;
   indots      = inpolygon(xdots,ydots,xrub,yrub);
   xindots     = xdots(find(indots > 0));
   yindots     = ydots(find(indots > 0));
   jump        = 0;
   for j = 1 : length(xindots)-1,
      if abs((abs(a_cuts) < 1) * (xindots(j) - xindots(j+1)) + ...
            (abs(a_cuts) >= 1) * (yindots(j) - yindots(j+1))) > 0.11,
         jump = [jump, j];
      end;
   end;
   if length(jump) > 1,
      if jump(2) == 1, jump(2) = []; end;
   end;
   jump = [jump, length(xindots)+1];
   diameter(i)   = 0;
   for j = 1 : length(jump)-1,
      diameter(i) = diameter(i) + sqrt((xindots(jump(j)+1) - xindots(jump(j+1)-1))^2 +...
         (yindots(jump(j)+1) - yindots(jump(j+1)-1))^2);
   end;
   
	%===================================================================================
   %  AREA_CROSS
   
   radius = (diameter(i)^2 / 4 + height(i)^2) / (2 * height(i));
   alpha  = 2 * atan2(diameter(i)/2, radius - height(i));
   area_cross(i) = 0.5 * radius^2 * (alpha - sin(alpha));
   
   %===================================================================================
   %  PRESTRESS
   
   prestress(i) = sumforce(i) / area_cross(i);
   
   
   else
   prestress(i) = NaN;
   end; %(if ~isempty(cutleft))
end; %(for i)


%  FIGURES

   figure(3);
   rect = [100, 100, 800, 560]; 
    load('MyColormaps','mycmap')
      set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;
%      set(gcf,'color',[1 1 1],'Position',rect,'Renderer','zbuffer'); clf; hold on;
     set(figure(3),'Colormap',mycmap);
%    set(gca,'color',[0 0 0]);
  

   caxis([0 1000])
   surf(x,y,sqrt(tx.^2+ty.^2)*1e12);
   whitebg('black');

%     view(2); colormap jet; shading interp; 
     view(2); shading interp;
     cbh = colorbar;
   set(gca,'FontSize',18,'FontWeight','bold'); box on;
    set(cbh,'FontSize',18,'FontWeight','bold');
   m3 = max2(sqrt(tx.^2+ty.^2)*1e12);
   h = plot3(xrub,yrub,m3*ones(size(xrub)),'w-');
   set(h,'LineWidth',2.0);
   h = quiver3(x,y,m3*ones(size(x)),tx,ty,zeros(size(x)),1);
   set(h,'Color','w','LineWidth',1);
      axis image; 
      axis tight; 
      axis ij;
%     axis([-115  115  -115  115]);
    axis([min(xvorig) max(xvorig) min(yvorig) max(yvorig)]);
%   set (gca, 'XTickLabel',{'','','','',''},'YTickLabel',{'','','','',''});
  xlabel('x (\mum)','FontSize',18,'FontWeight','bold'); ylabel('y (\mum)','FontSize',18,'FontWeight','bold');
tractime = str2num(savefilename(14:16))*10;
    headertrac = strcat(tractitle,';',savefilename(1:6),';',...
        num2str(tractime),'min');
   title(headertrac,'FontSize',18,'FontWeight','bold');   
   
%  Save the constrained traction figure to a frame 

   aviobj1 = addframe(aviobj1,figure(3));  

   figure(1);
    h = plot3(xrub,yrub,m1*ones(size(xrub)),'w-'); set(h,'LineWidth',2.0);

   aviobj3 = addframe(aviobj3,figure(1));  
           
   figure(2);
    h = plot3(xrub,yrub,m2*ones(size(xrub)),'w-'); set(h,'LineWidth',2.0);
   aviobj2 = addframe(aviobj2,figure(2));
   
%    figure(4);
%    left = 0.5;
%    h = text(left+mov,top+difv,'Constrained FFTC'); set(h,'FontSize',18,'FontWeight','bold');
%    text(left,top-difv,'RMS traction (Pa):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-2*difv,num2str(rmst_iterative*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-3*difv,'Orientation of principle tractions (^o):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-4*difv,num2str(theta0*180/pi,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-5*difv,'Net contractile moment (pNm):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-6*difv,num2str(-Trace_moment*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-7*difv,'Total strain energy (pJ):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-8*difv,num2str(Uecm*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-9*difv,'Max cumulative force (nN):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-10*difv,num2str(max(sumforce)*1e-3,'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-11*difv,'Prestress (Pa):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-12*difv,num2str(stnanmean(prestress),'%10.5f'),'FontSize',18,'FontWeight','bold');
%    text(left,top-13*difv,'Projected area of the cell (\mum^2):','FontSize',18,'FontWeight','bold'); 
%    text(left+mov,top-14*difv,num2str(area_cell,'%10.5f'),'FontSize',18,'FontWeight','bold');
  
 
%====================================================================================
%  Save program outputs

Pos      = str2num(savefilename(4:6));
Side     = str2num(savefilename(9:11));
time     = str2num(savefilename(14:15));

cellname = num2str(cellno);
%key = str2num(strcat(savefilename(4:6),savefilename(9:11),savefilename(14:15),num2str(cellname)));

% Entries in traction file are: 1. RMS traction, 2. Median traction, 3. CM, 4. Total force, and 5. Cell area      
tdata          = [Pos, Side, time, cellno, rmst_iterative*1e12,median_trac*1e12,-Trace_moment*1e12,tot_force*1e9,area_cell,-lambda(1,1)*1e12,-lambda(2,2)*1e12,theta*180/pi];

eval(['cd C:\Singlecell\traction\' cellfolder]);                                           
dlmwrite(fName,tdata,'-append',...  %# Print the matrix
         'delimiter','\t',...
         'newline','pc');
   
%====================================================================================

%cd('/Singlecell/Traction');
%filenametrac=strcat('/HTSBart/Traction/',cellfolder ,'/Traction_',savefilename(1:6),'_',num2str(cellno),'.dat');
%save(filenametrac,'tdata','-ASCII');
           
           end; %(cellno)
    end;
      
 aviobj1 = close(aviobj1);
 aviobj2 = close(aviobj2);
 aviobj3 = close(aviobj3);


clear a* A* B* c* C* d* D* h* i* j* k* m* n* p* r* s* T* t* u* U* v* x* y*;
close all;

return;