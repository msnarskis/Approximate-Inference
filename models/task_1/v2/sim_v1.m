%% Simulating for max_diff between conditions
clc; clear all;

% stimulus parameters
stim.eps_range = [0, 10];
stim.n = 25;
stim.ks = [3,5,7];
stim.kslab = {'k=3','k=5','k=7','k=9'};

stim.s = size(stim.ks,2);

% model parameters
par.var_t = 10; % (std dev)
par.var_n = 10;
par.var_v = 1;
par.var_sa = 100;
par.var_sv = 100;

par.pr_R = 0.5;
par.pr_C = 0.5;

par.nsamp = Inf; % trials per W vector
par.ntrials = 10000;

par.noisy_in = 1;

% Simulate and Analysis

resps = []; % for storing results

time = []; % comp time

for k=stim.ks
    sprintf("k=%d",k) % UX GUI
    sTime = cputime;

    stim.k = k;
    resp.k = k;

    stim = gen_stim_v1(stim);

    resp.center = model_v1(stim,par,0);
    resp.match  = model_v1(stim,par,1);

    time = [time, (cputime - sTime)] % sim done timestamp

    % correct
    resp.center_corr = abs(stim.corr - resp.center);
    resp.match_corr  = abs(stim.corr - resp.match);

    % by R and mean
    resp.match_corr_r1 = resp.match_corr(:,1);
    resp.match_corr_r0 = resp.match_corr(:,2);
    resp.match_r1 = resp.match(:,1);
    resp.match_r0 = resp.match(:,2);

    resp.center_corr_r1 = resp.center_corr(:,1);
    resp.center_corr_r0 = resp.center_corr(:,2);
    resp.center_r1 = resp.center(:,1);
    resp.center_r0 = resp.center(:,2);

    resp.match_corr_mean = mean([resp.match_corr_r1, resp.match_corr_r0],2);
    resp.center_corr_mean = mean([resp.center_corr_r1, resp.center_corr_r0],2);
    resp.match_mean = mean([resp.match_r1, resp.match_r0],2);
    resp.center_mean = mean([resp.center_r1, resp.center_r0],2);

    % diff by k
    resp.diff_corr = resp.match_corr-resp.center_corr;
    resp.diff_corr_mean = resp.match_corr_mean - resp.center_corr_mean;
    resp.diff = resp.match-resp.center;
    resp.diff_mean = resp.match_mean - resp.center_mean;

    resps = [resps; resp];
end

%% Save
equiprob = 0;
if equiprob
    file = strcat(pwd, sprintf("/models/task_1/sims/simS-max-diff-class_t-%d_k-%d_nsamp-%d_pC%1.1f.mat", par.ntrials, stim.ks(end), par.nsamp, par.pr_C));
else
    file = strcat(pwd, sprintf("/models/task_1/v2/sims/sim-max-diff-class_t-%d_k-%d_nsamp-%d_pC%1.1f.mat", par.ntrials, stim.ks(end), par.nsamp, par.pr_C));
end
save(file,'resps','par','stim');

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
g_end = 0.7;
b_start= 1;
b_end = 0;
K_COL_bg = [linspace(r_start, r_end, s); linspace(g_start, g_end, s); linspace(b_start, b_end, s); linspace(0.6,0.6,s)]';

% misc
position_huge = [100 0 1700 900];
position_ok = [100 200 1600 600];
position_ppt = [100 100 900 450;];

% plot psychometric curves and diff between match, center
figure
hold on

% for legends
lbl1 = [];
lbl2 = [];
lbl3 = [];

% loop over k's
for i=1:s
    % match
    subplot(s,3,3*i-2);hold on;
    plot(stim.eps, resps(i).match_corr, '--', 'Color', K_COL_BG(i,:), 'LineWidth', 1);
    lbl1 =plot(stim.eps, resps(i).match_corr_mean, 'Color', K_COL_BG(i,:),...
			'LineWidth', 3, 'DisplayName', sprintf("matched", stim.ks(i)));
		
	% center
    plot(stim.eps, resps(i).center_corr, '--', 'Color', K_COL_bg(i,:), 'LineWidth', 1);
    lbl2 = plot(stim.eps, resps(i).center_corr_mean, 'Color', K_COL_bg(i,:),...
			'LineWidth', 3, 'DisplayName', sprintf("center", stim.ks(i)));
	
	title(sprintf("k=%d", stim.ks(i)));
	xlim(stim.eps_range);ylim([0.5,1]);
	xlabel("stim: dist from center");ylabel("P(correct)");
	legend([lbl1 lbl2], 'Location','best');

    % difference
    subplot(s,3,[2 3*s]);hold on;
    plot(stim.eps, resps(i).diff_corr, '--', 'Color', K_COL_BG(i,:), 'LineWidth', 1);
    lbl3 = [lbl3 plot(stim.eps, resps(i).diff_corr_mean, 'Color', K_COL_BG(i,:),...
			'LineWidth', 5, 'DisplayName', sprintf("k=%d", stim.ks(i)))];
end

subplot(s,3,1);
title(sprintf("Majority Task Psychometric Curve, %d sample -- k=%d",par.nsamp,stim.ks(i)));

% subplot(2,2,3);
% xlim(stim.eps_range);
% title(sprintf("Psychometric: Center -- nsamp:%d",par.nsamp));
% xlabel("stim: dist from center");ylabel("P(correct)");
% legend(lbl2,stim.kslab,'Location','best');

subplot(s,3,[2 3*s]);ylim([-.01,.06]);
xlim(stim.eps_range);
title(sprintf("Difference in curves (M-C) -- %d sample",par.nsamp));
xlabel("stim: dist from center");ylabel("\Delta P(correct)");
legend(lbl3,stim.kslab,'Location','best');

set(gcf,'Position',position_ppt);


%% plot flattened curves max diff
% 
% for i=1:s
%     f = figure;
%     clim = [0 1];
%     xlims = [0.5, size(resps(i).match_corr,1)+0.5];
%     ylims = [0.5, size(resps(i).match_corr,2)+0.5];
% 
%     % match condition
%     ax1 = subplot(2,2,1);
%     hold on
% 
%     imagesc(resps(i).match_corr', clim);
%     xlabel("stimulus location");ylabel("stimulus class");
%     xlim(xlims);ylim(ylims);yticks([]);
%     title(sprintf("Flattened psychometric curves: Match; k=%d, nsamp=%d",stim.ks(i),par.nsamp));
%     c = colorbar; ylabel(c, "P(correct)");
% 
%     % center condition
%     ax2 = subplot(2,2,3);
%     hold on
% 
%     imagesc(resps(i).center_corr', clim);
%     xlabel("stimulus location");ylabel("stimulus class");
%     xlim(xlims);ylim(ylims);yticks([]);
%     title(sprintf("Flattened psychometric curves: Center; k=%d, nsamp=%d",stim.ks(i),par.nsamp));
%     c = colorbar; ylabel(c, "P(correct)");
% 
%     % difference between condsd
%     ax3 = subplot(2,2,[2,4]);
%     hold on
% 
%     imagesc(resps(i).diff_corr');
%     xlabel("stimulus location");ylabel("stimulus class");
%     xlim(xlims);ylim(ylims);yticks([]);
%     title(sprintf("Difference btw Match - Center; k=%d, nsamp=%d",stim.ks(i),par.nsamp));
%     c = colorbar; ylabel(c, "\Delta P(correct)");
% 
%     set(gcf,'Position',position_huge);
% end
% 
%% Diff scaling by k

diff_mean = zeros(1,s);
diff_part = nan(1,s);

for i=1:s
    diff_mean(i) = max(resps(i).diff_corr_mean);
    diff_part(i) = max(max(resps(i).diff_corr));
end

figure;
hold on

p1 = plot(stim.ks,diff_mean, 'LineWidth', 4);
p2 = plot(stim.ks,diff_part, 'LineWidth', 4);

legend([p1, p2], {'Average', 'By Stimulus Class'}, 'Location', 'northwest');
title(sprintf("Maximum Diff in P(Correct) by k, nsamp=%d",par.nsamp));
ylabel("Max \Delta P(correct)");xlabel("k");
xticks(stim.ks)

%% PLot more
if 0
load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_1\v2\sims\sim-max-diff-class_t-10000_k-7_nsamp-Inf_pC0.5.mat')
resps1=resps;
par1 = par;
stim1 = stim;

load('C:\Users\User\MATLAB\Projects\R1-Approx-Inf\models\task_1\v2\sims\sim-max-diff-class_t-10000_k-7_nsamp-1_pC0.5.mat')
respsInf = resps;
parInf = par;
stimInf = stim;

for i=1:s
    lbl1 = [];
    lbl2 = [];
    lbl3 = [];
    lbl4 = [];
    figure;
    % match & center 1 samp
    subplot(2,2,1);hold on;
    plot(stim1.eps, resps1(i).match_corr, '--', 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 1);
    lbl1 = [lbl1 plot(stim1.eps, resps1(i).match_corr_mean, 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 5)...
        plot(stim1.eps, resps1(i).center_corr_mean, 'Color',[i/s, 1-(i/s), 0.7,], 'LineWidth', 5)];

    xlim(stim1.eps_range);
    title(sprintf("Psychometric: k=%d nsamp:%d",stim1.ks(s),par1.nsamp));
    xlabel("stim: dist from center");ylabel("P(correct)");
    legend(lbl1,stim1.kslab,'Location','best');

    % match & center Inf samp
    subplot(2,2,2);hold on;
    plot(stimInf.eps, respsInf(i).match_corr, '--', 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 1);
    lbl2 = [lbl2 plot(stimInf.eps, respsInf(i).match_corr_mean, 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 5)...
        plot(stimInf.eps, respsInf(i).center_corr_mean, 'Color',[i/s, 1-(i/s), 0.7], 'LineWidth', 5)];

    xlim(stimInf.eps_range);
    title(sprintf("Psychometric: k=%d nsamp:%d",stimInf.ks(s),parInf.nsamp));
    xlabel("stim: dist from center");ylabel("P(correct)");
    legend(lbl2,stimInf.kslab,'Location','best');

    % difference 1 samp
    subplot(2,2,3);hold on;
    plot(stim1.eps, resps1(i).diff_corr, '--', 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 1);
    lbl3 = [lbl3 plot(stim1.eps, resps1(i).diff_corr_mean, 'Color',[i/s,0.3, 1-(i/s)], 'LineWidth', 5)];

    xlim(stim1.eps_range);
    title(sprintf("Difference in curves (M-C) -- k=%d nsamp:%d",stim1.ks(s), par1.nsamp));
    xlabel("stim: dist from center");ylabel("\Delta P(correct)");
    legend(lbl3,stim1.kslab,'Location','best');

    % difference Inf samp
    subplot(2,2,4);hold on;
    plot(stimInf.eps, respsInf(i).diff_corr, '--', 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 1);
    lbl4 = [lbl4 plot(stimInf.eps, resps(i).diff_corr_mean, 'Color',[i/s,0.3, 1-(i/s)], 'LineWidth', 5)];

    xlim(stimInf.eps_range);
    title(sprintf("Difference in curves (M-C) -- k=%d nsamp:%d",stimInf.ks(s), parInf.nsamp));
    xlabel("stim: dist from center");ylabel("\Delta P(correct)");
    legend(lbl4,stimInf.kslab,'Location','best');

    set(gcf,'Position',[100 100 1600 1000]);
end
end

%% Aux Functions

% takes a properly weighted average of stim classes
% for full experiment
function [exp] = class_collapse(resp,k,w)
    if w
        % init
        exp = zeros(size(resp,1),1);

        % loop over stim classes
        for i=ceil(k/2)
            exp = exp + resp(:,i)*nchoosek(k,i-1);
        end

        % norm
        exp = exp/2^k;
    else
        exp = mean(resp,2);
    end
end
