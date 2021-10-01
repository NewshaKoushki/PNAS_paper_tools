eval(['cd C:\Singlecell\traction']); pwd; 

load('all_scatter_mutant_26.mat')
load('all_scatter_ctrl_26.mat')

load('MyRedColormap','mycmap_red')
load('MyBlueColormap','mycmap_blue')

T_Ave_v_mutant_26 = cat(1,T_Ave_v_103,T_Ave_v_111,T_Ave_v_3113,T_Ave_v_3118,T_Ave_v_3120);
H_Cum_v_mutant_26 = cat(1,H_Cum_v_103,H_Cum_v_111,H_Cum_v_3113,H_Cum_v_3118,H_Cum_v_3120);

T_Ave_v_ctrl_26 = cat(1,T_Ave_v_104,T_Ave_v_112,T_Ave_v_3123,T_Ave_v_3125,T_Ave_v_3129);
H_Cum_v_ctrl_26 = cat(1,H_Cum_v_104,H_Cum_v_112,H_Cum_v_3123,H_Cum_v_3125,H_Cum_v_3129);

%% scatter plot

figure(1),

% plot the ctrl in red, different symbols for different cells
scatter(T_Ave_v_103,H_Cum_v_103,'r','.');axis([0 4000 0 150]), hold on, 
scatter(T_Ave_v_111,H_Cum_v_111,'r','+'); 
scatter(T_Ave_v_3113,H_Cum_v_3113,'r','x'); 
scatter(T_Ave_v_3118,H_Cum_v_3118,'r','*');
scatter(T_Ave_v_3120,H_Cum_v_3120,'r','o');

% plot the ctrl in blue, different symbols for different cells
scatter(T_Ave_v_104,H_Cum_v_104,'b','.'); 
scatter(T_Ave_v_112,H_Cum_v_112,'b','+'); 
scatter(T_Ave_v_3123,H_Cum_v_3123,'b','x'); 
scatter(T_Ave_v_3125,H_Cum_v_3125,'b','*');
scatter(T_Ave_v_3129,H_Cum_v_3129,'b','o');

%% 2d histogram

bin_time = (5:10:145);
bin_force = (100:200:3900);

data_ctrl_26 = cat(2,T_Ave_v_ctrl_26,H_Cum_v_ctrl_26);
data_mutant_26 = cat(2,T_Ave_v_mutant_26,H_Cum_v_mutant_26);

n_ctrl_26 = hist3(data_ctrl_26, {bin_force,bin_time});
nH1 = n_ctrl_26';
nH1(size(n_ctrl_26,2)+1,size(n_ctrl_26,1)+1) = 0;
yb = linspace(0,150,size(n_ctrl_26,2)+1);
xb = linspace(0,4000,size(n_ctrl_26,1)+1);

figure (2),
h = contourf(xb,yb,nH1/sum(sum(nH1))*100,150,'EdgeColor', 'none'); 
colorbar
axis([0 3000 0 150])
colormap (mycmap_blue)
view(2)
box off 
grid off
caxis ([0 5]);
xlabel('Average Traction');
ylabel('Persistence');


n_mutant_26 = hist3(data_mutant_26, {bin_force,bin_time});
nH1 = n_mutant_26';
nH1(size(n_mutant_26,2)+1,size(n_mutant_26,1)+1) = 0;
yb = linspace(0,150,size(n_mutant_26,2)+1);
xb = linspace(0,4000,size(n_mutant_26,1)+1);

figure (3),
h = contourf(xb,yb,nH1/sum(sum(nH1))*100,150,'EdgeColor', 'none'); 
colorbar
axis([0 3000 0 150])
colormap (mycmap_red)
view(2)
box off 
grid off
caxis ([0 5]);
xlabel('Average Traction (Pa)');
ylabel('Persistence time (hours)');


