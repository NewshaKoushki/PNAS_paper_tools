
pos = 21;
time_interval = 22; % interval between frames in mins
% Get result file names
[nameMat, pathnameMat] = uigetfile({'*.mat','Matlab data file (*.mat)'}, 'Select Mat file for YAP intensity');
if ~nameMat, warning('No file selected.'); return; end
filenameYAP = [pathnameMat nameMat];

[nameTraction, pathnameTraction] = uigetfile({'*.mat','Matlab data file (*.mat)'}, 'Select the Mat file for Traction force results');
if ~nameTraction, warning('No file selected.'); return; end
filenameTraction = [pathnameTraction nameTraction];


% Load Mat files
load(filenameYAP);
load(filenameTraction);

%% Graph2 Nucleus YAP & Strain Energy vs time
fig1 = figure;

for i = 1:size(NucYAP(pos,:),2)
    
    % plot on the left
    plot([0:i-1]*time_interval,NucYAP(pos,1:i),'bo','MarkerFaceColor','b','markers',3);
    h1 = gca;
    set(h1, 'YAxisLocation', 'left', 'xlim', get(h1, 'xlim'));
    xlabel('time (min)','fontsize',16);
    ylabel( 'YAP intensity nucleus','fontsize',16); 
    xlim([0 (size(NucYAP(pos,:),2)-1)*time_interval]);
    ylim([0 100000]);
    h1 = gca;
    
    % plot on the right
    h_ax_line = axes('position', get(h1, 'position')); % Create a new axes in the same position as the first one, overlaid on top
    set(h_ax_line, 'YAxisLocation', 'right', 'xlim', get(h1, 'xlim'), 'color', 'none'); % Put the new axes' y labels on the right, set the x limits the same as the original axes', and make the background transparent
    
 hold
    plot(h_ax_line,[0:i-1]*time_interval,StrainEnergy(pos,1:i),'ro','MarkerFaceColor','r','markers',3);
   
  
    ylabel(h_ax_line, 'Strain Energy (p)','fontsize',16);
    ylim(h_ax_line,[0 2]);
  
    pause (0.5)
    frame = getframe(1);
    util.write_gif(util.gen_addr([pathnameMat '\' 'pos' num2str(pos) '_NucYAPvsStrainEnergy.gif']), 0.5 + (i == 1)*1.5, frame.cdata);
    
end

frame2 = getframe(1);
%export_fig(frame2,util.gen_addr([pathnameMat '\' 'pos' num2str(pos) '_Plot_NucYAPvsStrain' '.tif']),'-transparent');

%% Graph2 Nucleus YAP & rmst vs time

fig2 = figure;

for i = 1:size(NucYAP(pos,:),2)
    
    % plot on the left
    plot([0:i-1]*time_interval,NucYAP(pos,1:i),'bo','MarkerFaceColor','b','markers',3);
    h1 = gca;
    set(h1, 'YAxisLocation', 'left', 'xlim', get(h1, 'xlim'));
    xlabel('time (min)','fontsize',16);
    ylabel( 'YAP intensity nucleus','fontsize',16); 
    xlim([0 (size(NucYAP(pos,:),2)-1)*time_interval]);
    ylim([0 100000]);
    h1 = gca;
    
    % plot on the right
    h_ax_line = axes('position', get(h1, 'position')); % Create a new axes in the same position as the first one, overlaid on top
    set(h_ax_line, 'YAxisLocation', 'right', 'xlim', get(h1, 'xlim'), 'color', 'none'); % Put the new axes' y labels on the right, set the x limits the same as the original axes', and make the background transparent
    
 hold
    plot(h_ax_line,[0:i-1]*time_interval,rmst(pos,1:i),'ro','MarkerFaceColor','r','markers',3);
   
  
    ylabel(h_ax_line, 'RMST (Pa)', 'fontsize',16);
    ylim(h_ax_line,[0 1200]);
  
    pause (0.5)
    frame = getframe(1);
    util.write_gif(util.gen_addr([pathnameMat '\' 'pos' num2str(pos) '_NucYAPvsRMST.gif']), 0.5 + (i == 1)*1.5, frame.cdata);
    
end

 frame1 = getframe(1);
%export_fig(frame1,util.gen_addr([pathnameMat '\' 'pos' num2str(pos) '_Plot_NucYAPvsRMST' '.tif']),'-transparent');


