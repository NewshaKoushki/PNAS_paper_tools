eval(['cd C:\Singlecell\traction']); pwd; 

load('all_scatter_mutant_26.mat')
load('all_scatter_ctrl_26.mat')

load('MyRedColormap','mycmap_red')
load('MyBlueColormap','mycmap_blue')
load('MyColormaps','mycmap')
load('MyBluetoRedColormap','mycmap_bluetored')

TH_Cum_v_mutant_26 = cat(1,TH_Cum_v_103,TH_Cum_v_111,TH_Cum_v_3113,TH_Cum_v_3118,TH_Cum_v_3120);
H_Cum_v_mutant_26 = cat(1,H_Cum_v_103,H_Cum_v_111,H_Cum_v_3113,H_Cum_v_3118,H_Cum_v_3120);


TH_Cum_v_ctrl_26 = cat(1,TH_Cum_v_104,TH_Cum_v_112,TH_Cum_v_3123,TH_Cum_v_3125,TH_Cum_v_3129);
H_Cum_v_ctrl_26 = cat(1,H_Cum_v_104,H_Cum_v_112,H_Cum_v_3123,H_Cum_v_3125,H_Cum_v_3129);




TH_Ave_v_mutant_26 = (TH_Cum_v_mutant_26 ./ H_Cum_v_mutant_26);
TH_Ave_v_ctrl_26 = (TH_Cum_v_ctrl_26 ./ H_Cum_v_ctrl_26);


    

%% scatter plot

% figure(1),
% 
% % plot the ctrl in red, different symbols for different cells
% scatter(T_Cum_v_103,H_Cum_v_103,'r','.');axis([0 4000 0 150]), hold on, 
% scatter(T_Cum_v_111,H_Cum_v_111,'r','+'); 
% scatter(T_Cum_v_3113,H_Cum_v_3113,'r','x'); 
% scatter(T_Cum_v_3118,H_Cum_v_3118,'r','*');
% scatter(T_Cum_v_3120,H_Cum_v_3120,'r','o');
% 
% % plot the ctrl in blue, different symbols for different cells
% scatter(T_Cum_v_104,H_Cum_v_104,'b','.'); 
% scatter(T_Cum_v_112,H_Cum_v_112,'b','+'); 
% scatter(T_Cum_v_3123,H_Cum_v_3123,'b','x'); 
% scatter(T_Cum_v_3125,H_Cum_v_3125,'b','*');
% scatter(T_Cum_v_3129,H_Cum_v_3129,'b','o');

figure(2),
scatter((TH_Ave_v_mutant_26),(H_Cum_v_mutant_26),'r','.'), hold on
scatter((TH_Ave_v_ctrl_26),(H_Cum_v_ctrl_26),'b','.');
axis([0 6000 0 150])
%% 2d histogram

bin_time = (2.5:5:142.5);
bin_force = (100:200:5900);

data_ctrl_26 = cat(2,(TH_Ave_v_ctrl_26),H_Cum_v_ctrl_26);
data_mutant_26 = cat(2,(TH_Ave_v_mutant_26),H_Cum_v_mutant_26);




n_ctrl_26 = hist3(data_ctrl_26, {bin_force,bin_time});
nH1 = n_ctrl_26';
nH1(size(n_ctrl_26,2)+1,size(n_ctrl_26,1)+1) = 0;
yb = linspace(0,150,size(n_ctrl_26,2)+1);
xb = linspace(0,6000,size(n_ctrl_26,1)+1);



figure(3),
h1 = contourf(xb,yb,nH1*100/sum(sum((nH1))),150,'EdgeColor', 'none'); 
colorbar
%axis([7 13 0 150])
colormap (mycmap_blue)
view(2)
box off 
grid off
caxis ([0 6]);
xlabel('Average Traction (Pa)');
ylabel('Persistence time (hours)');
plot1 = frame2im(getframe);



n_mutant_26 = hist3(data_mutant_26, {bin_force,bin_time});
nH2 = n_mutant_26';
nH2(size(n_mutant_26,2)+1,size(n_mutant_26,1)+1) = 0;
yb = linspace(0,150,size(n_mutant_26,2)+1);
xb = linspace(0,6000,size(n_mutant_26,1)+1);

figure(4),
h2 = contourf(xb,yb,nH2*100/sum(sum((nH2))),150,'EdgeColor', 'none');  
colorbar
%axis([7 13 0 150])
colormap(mycmap_red)
view(2)
box off 
grid off
caxis ([0 6]);
xlabel('Average Traction (Pa)');
ylabel('Persistence time (hours)');
plot2 = frame2im(getframe);

figure(5)
image(plot1);
alpha(0.3);
hold on
image(plot2);
alpha(0.3);

%%

figure(6),
h3 = contourf(xb,yb,nH2*100/sum(sum((nH2)))-nH1*100/sum(sum((nH1))),150,'EdgeColor', 'none'); 
colorbar
%axis([7 13 0 150])
colormap (mycmap_bluetored)
view(2)
box off 
grid off
caxis ([-5 5]);
xlabel('Average Traction (Pa)');
ylabel('Persistence time (hours)');

%% get average values for hotspot persistence and traction in hotspots

%ctrl
HS_life_ave_ctrl = 0;
THS_ave_ctrl = 0;
HS_number_ctrl = length(H_Cum_v_ctrl_26)-sum(isnan(H_Cum_v_ctrl_26));

for i=1:length(H_Cum_v_ctrl_26)
    if isnan(H_Cum_v_ctrl_26(i)) == 0,
    HS_life_ave_ctrl = HS_life_ave_ctrl + H_Cum_v_ctrl_26(i);
    THS_ave_ctrl = THS_ave_ctrl + TH_Ave_v_ctrl_26(i);
    end
    
end

HS_life_ave_ctrl = HS_life_ave_ctrl/HS_number_ctrl;
THS_ave_ctrl = THS_ave_ctrl/HS_number_ctrl;

%mutants
HS_life_ave_mutant = 0;
THS_ave_mutant = 0;
HS_number_mutant = length(H_Cum_v_mutant_26)-sum(isnan(H_Cum_v_mutant_26));

for i=1:length(H_Cum_v_mutant_26)
    if isnan(H_Cum_v_mutant_26(i)) == 0,
    HS_life_ave_mutant = HS_life_ave_mutant + H_Cum_v_mutant_26(i);
    THS_ave_mutant = THS_ave_mutant + TH_Ave_v_mutant_26(i);
    end
    
end

HS_life_ave_mutant = HS_life_ave_mutant/HS_number_mutant;
THS_ave_mutant = THS_ave_mutant/HS_number_mutant;

