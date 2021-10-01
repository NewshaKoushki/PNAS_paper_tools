% This code loads the .mat files that gather the 2d histogram info for each cell of
% the different groups, and plot one 2d histogram per group.

eval(['cd C:\Singlecell\traction']); pwd; 

load('all_histo_ctrl_4.mat')
load('all_histo_ctrl_26.mat')
load('all_histo_mutants_4.mat')
load('all_histo_mutants_26.mat')

%==========================================================================
% ctrl 26

histo_ctrl_26 = nH1_104 + nH1_112 + nH1_3123 + nH1_3125 + nH1_3129;

figure,
h = surf(xb,yb,-log((histo_ctrl_26/sum(sum(histo_ctrl_26)))),'EdgeColor', 'none'); 
colorbar
axis([0 145 0 7000])
colormap (hot)
view(2)
box off 
grid off
caxis ([4 11]);

%==========================================================================
% ctrl 4
histo_ctrl_4 = nH1_105 + nH1_3132 + nH1_3134 + nH1_3136 + nH1_3137;

figure,
h = surf(xb,yb,-log((histo_ctrl_4/sum(sum(histo_ctrl_4)))),'EdgeColor', 'none'); 
colorbar
axis([0 145 0 7000])
colormap (hot)
view(2)
box off 
grid off
caxis ([4 11]);

%==========================================================================
% mutant 26

histo_mutant_26 = nH1_103 + nH1_111 + nH1_3113 + nH1_3117 + nH1_3118 + nH1_3120;

figure,
h = surf(xb,yb,-log((histo_mutant_26/sum(sum(histo_mutant_26)))),'EdgeColor', 'none'); 
colorbar
axis([0 145 0 7000])
colormap (hot)
view(2)
box off 
grid off
caxis ([4 11]);

%==========================================================================
% mutant 4

histo_mutant_4 = nH1_102 + nH1_3101 + nH1_3102 + nH1_3103 + nH1_3110;

figure,
h = surf(xb,yb,-log((histo_mutant_4/sum(sum(histo_mutant_4)))),'EdgeColor', 'none'); 
colorbar
axis([0 145 0 7000])
colormap (hot)
view(2)
box off 
grid off
caxis ([4 11]);


