% FTTC calculates the tractions, using both the Unconstrained
%	and the Constrained FTTC.
%	The required input is a file of displacements (4 columns: x-coordinates,
%	y-coordinates, x-displacements, y-displacements, all in pixels), 
%	pixel to micron conversion factor, Young's modulus and Poisson's ratio 
%	of the gel, and the cell boundary file (it is needed for the 
%	Constrained FTTC only).

%	Iva Marija Tolic-Norrelykke 03-22-01
   
%====================================================================================
%  INPUT: LOAD DISPLACEMENTS

   clear all;
   %close all;
   [filename,pathname] = uigetfile('*.*','Input File: Displacements');
   savedirectory = pathname;
   savefilename  = filename;
   cd(pathname);
   displ = load(filename);

%    x_orig = displ(:,1);
%    y_orig = displ(:,3).^2 + displ(:,4).^2;
% 	
%    figure(10);hold on;
%    set(gca,'Color',[1.000 1.000 1.000]); 
%    title(strcat('Displacement'),'FontSize',20,'FontWeight','bold');grid on;
%    plot(x_orig,y_orig,'r.');
%    xlabel('Coord 1','FontSize',12,'FontWeight','bold'); ylabel('Disp','FontSize',12,'FontWeight','bold');
%    set(gca,'FontSize',12,'FontWeight','bold');
%    hold off; box on; drawnow;
%====================================================================================

%    return;

max_x = max(displ(:,1));
min_x = min(displ(:,1));

x_orig = displ(:,1);
y_orig = displ(:,3).^2 + displ(:,4).^2;

step_num = 1600;
step_x = (max_x-min_x)/step_num;
dispX_matrix = zeros(step_num+1,1);
dispY_matrix = ones(step_num+1,1).*-1;

% for k= 1 : step_num+1,
%      dispX_matrix(k) = min_x+(k-1)*step_x;
% end;

for k= 1 : size(x_orig,1),
     col_num = min(round((x_orig(k)-min_x)/step_x)+1,step_num+1);                                                                                                                                                 
     dispX_matrix(col_num)  =  x_orig(k);
     if dispY_matrix(col_num) < y_orig(k),
         dispY_matrix(col_num) = y_orig(k);
     end;
end;

%get rid of the zero elements;
disp_matrix = [];
for k= 1 : size(dispY_matrix,1),
     if dispY_matrix(k) > 0,
         disp_matrix = [disp_matrix; [dispX_matrix(k) dispY_matrix(k)] ];
     end;
end;

% dispX_matrix = [];
% dispY_matrix = [];
% for k= 1 : step_num+1,
%      dispX_matrix(k) = min_x+(k-1)*step_x;
% end;
% dispY_matrix = interp1(disp_matrix(:,1),disp_matrix(:,2),dispX_matrix(:),'spline');

model_x = disp_matrix(:,1)./10;
model_y = disp_matrix(:,2);
model_y = model_y./max(model_y);

figure(10);hold on;
set(gca,'Color',[1.000 1.000 1.000]); 
title(strcat('Displacement'),'FontSize',20,'FontWeight','bold');grid on;
plot(model_x,model_y,'-r');
axis([-10 10 0 +1.5]);
xlabel('Coord 1','FontSize',12,'FontWeight','bold'); ylabel('Disp','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold');
hold off; box on; drawnow;


save_file = strcat('G:\Research\AndrewBead3D\JimBulter1\Model_disp','.mat');
save (save_file, 'model_x', 'model_y');
	  

return;

