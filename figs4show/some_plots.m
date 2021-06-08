%% Load Data
% task 1
load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_1\v2\sims\sim-max-diff-class_t-10000_k-7_nsamp-1_pC0.5.mat')
t1resp1=resps;
t1par1 = par;
t1stim1 = stim;

load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_1\v2\sims\sim-max-diff-class_t-10000_k-7_nsamp-Inf_pC0.5.mat')
t1respInf=resps;
t1parInf = par;
t1stimInf = stim;

% task 2
load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_2\v2\sims\sim_max_diff_class_t-5000-k-7-nsamp-1-pR0-5.mat')
t2resp1=resps;
t2par1 = par;
t2stim1 = stim;

load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_2\v2\sims\sim_max_diff_class_t-5000-k-7-nsamp-Inf-pR0-5.mat')
t2respInf=resps;
t2parInf = par;
t2stimInf = stim;

% task 3
load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_3\sims\sim-max-diff-class_t-8000_k-7_nsamp-1_pC0.5.mat')
t3resp1=resps;
t3par1 = par;
t3stim1 = stim;

load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_3\sims\sim-max-diff-class_t-8000_k-7_nsamp-Inf_pC0.5.mat')
t3respInf=resps;
t3parInf = par;
t3stimInf = stim;

s = size(stim.ks,2);

%% calcs
t1comp_mean_1 = zeros(1,s-1);
t1comp_part_1 = zeros(1,s-1);
t1comp_mean_Inf = zeros(1,s-1);
t1comp_part_Inf = zeros(1,s-1);

t2comp_mean_1 = zeros(1,s);
t2comp_part_1 = zeros(1,s);
t2comp_mean_Inf = zeros(1,s);
t2comp_part_Inf = zeros(1,s);

t3comp_mean_1 = zeros(1,s);
t3comp_part_1 = zeros(1,s);
t3comp_mean_Inf = zeros(1,s);
t3comp_part_Inf = zeros(1,s);

for i=1:4
	if i < 4
		t1comp_mean_1(i) = max(t1resp1(i).diff_corr_mean);
		t1comp_part_1(i) = max(max(t1resp1(i).diff_corr));

		t1comp_mean_Inf(i) = max(t1respInf(i).diff_corr_mean);
		t1comp_part_Inf(i) = max(max(t1respInf(i).diff_corr));

	end

	% task 2
	t2comp_mean_1(i) = max(t2resp1(i).diff_corr_mean);
	t2comp_part_1(i) = max(max(t2resp1(i).diff_corr));
	
	t2comp_mean_Inf(i) = max(t2respInf(i).diff_corr_mean);
	t2comp_part_Inf(i) = max(max(t2respInf(i).diff_corr));
	
	% task 3
	t3comp_mean_1(i) = max(t3resp1(i).diff_corr_mean);
	t3comp_part_1(i) = max(max(t3resp1(i).diff_corr));
	
	t3comp_mean_Inf(i) = max(t3respInf(i).diff_corr_mean);
	t3comp_part_Inf(i) = max(max(t3respInf(i).diff_corr));
end

%% configure for plots
s = stim.s;
% colors - index by: K_COL(i,:)
r_start = 0;
r_end = 1;
g_start = 0.15; 
g_end = 0.3;
b_start= 1;
b_end = 0;
K_COL_BG = [linspace(r_start, r_end, s); linspace(g_start, g_end, s); linspace(b_start, b_end, s)]';

r_start = 0;
r_end = 1;
g_start = 0.6; 
g_end = 0;
b_start= 1;
b_end = 0;
K_COL_bg = [linspace(r_start, r_end, s); linspace(g_start, g_end, s); linspace(b_start, b_end, s); linspace(0.6,0.6,s)]';

% misc
position_huge = [100 0 1700 900];
position_ok = [100 200 1600 600];
position_ppt = [100 100 900 450;]


%% plot full model compounding effects
figure;
hold on

plot(stim.ks(1:3), t1comp_mean_1, 'Color', '#D95319', 'LineWidth', 3,'DisplayName', 'Majority (max avg)');
plot(stim.ks(1:3), t1comp_part_1, ':', 'Color', '#D95319', 'LineWidth', 3,'DisplayName', 'Majority (max by class)');

plot(stim.ks, t2comp_mean_1, 'Color','#0072BD', 'LineWidth', 3,'DisplayName', 'All-or-None (max avg)');
plot(stim.ks, t2comp_part_1, ':', 'Color', '#0072BD', 'LineWidth', 3,'DisplayName', 'All-or-None (max by class)');

plot(stim.ks, t3comp_mean_1, 'Color', 	'#EDB120', 'LineWidth', 3,'DisplayName', 'Even/Odd (max avg)');
plot(stim.ks, t3comp_part_1, ':', 'Color', '#EDB120', 'LineWidth', 3,'DisplayName', 'Even/Odd (max by class)');

title("Compounding of Improvement for Full Models");
xlabel("number of frames (k)");
ylabel("maximum improvement in performance");
legend('Location',"eastoutside")










