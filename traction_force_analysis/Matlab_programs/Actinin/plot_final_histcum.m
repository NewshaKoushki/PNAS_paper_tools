% This code loads the .mat files that gather the histogram info to represent 
% the hotspots persistency for each cell of the different groups, and plot 
% 2 histograms to compare the groups at different stiffnesses. 

eval(['cd C:\Singlecell\traction']); pwd; 

load('all_histocum_ctrl_4.mat')
load('all_histocum_ctrl_26.mat')
load('all_histocum_mutants_4.mat')
load('all_histocum_mutants_26.mat')

% ctrl 4
hist_cum_ctrl_4 = hist_cum_105 + hist_cum_3132 + hist_cum_3134 + hist_cum_3136 + hist_cum_3137;

% ctrl 26
hist_cum_ctrl_26 = hist_cum_104 + hist_cum_112 + hist_cum_3123 + hist_cum_3125 + hist_cum_3129;

%mutant 26
hist_cum_mutant_26 = hist_cum_103 + hist_cum_111 + hist_cum_3113 + hist_cum_3117 + hist_cum_3118 + hist_cum_3120;

%mutant 4
hist_cum_mutant_4 = hist_cum_3101 + hist_cum_3102 +hist_cum_3103 + hist_cum_102 + hist_cum_3110;

figure(1),
bwt=bar(bin_time,log(hist_cum_ctrl_26/sum(hist_cum_ctrl_26)));

chw = get(bwt,'child');
set(chw,'facea',.5,'FaceColor',[1 0 0],'EdgeColor','k');
 
hold on
bk=bar(bin_time,log(hist_cum_mutant_26/sum(hist_cum_mutant_26)));
chk = get(bk,'child');
set(chk,'facea',.5,'FaceColor',[0 0 1],'EdgeColor','k');
% axis([0 150 0 0.5])
% 
xlabel('Time(x10 minutes)');
ylabel('log(Count/Sum)');
% 
% 
%  

figure(2),
bwt=bar(bin_time,log(hist_cum_ctrl_4/sum(hist_cum_ctrl_4)));

chw = get(bwt,'child');
set(chw,'facea',.5,'FaceColor',[1 0 0],'EdgeColor','k');
 
hold on
bk=bar(bin_time,log(hist_cum_mutant_4/sum(hist_cum_mutant_4)));
chk = get(bk,'child');
set(chk,'facea',.5,'FaceColor',[0 0 1],'EdgeColor','k');
% axis([0 150 0 0.5])
% 
xlabel('Time(x10 minutes)');
ylabel('log(Count/Sum)');
