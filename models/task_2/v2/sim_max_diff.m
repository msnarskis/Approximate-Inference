%% Simulating for max_diff between conditions
clc, clear all;

% stimulus parameters
stim.eps_range = [0, 8];
stim.n = 30;
stim.ks = [2,3,5,7];
stim.kslab = {'k=2','k=3','k=5','k=7'};

% model parameters
par.var_t = 10; % (std dev)
par.var_n = 10;
par.var_v = 1;
par.var_sa = 100;
par.var_sv = 100;

par.pr_R = 0.5;
par.pr_C = 0.5;

par.nsamp = 1; % trials per W vector
par.ntrials = 100;

par.noisy_in = 0;

%% Simulate
model = Model_Task_2(par);

resps = [];

for k=stim.ks
    stim.k = k;
    resp.k = k;
    
    stim = generate_stim_v2_2(stim);
    
    resp.center = model.simulate(stim,0);
    resp.match  = model.simulate(stim,1);
    
    resp.center_corr = abs(stim.corr - resp.center);
    resp.match_corr  = abs(stim.corr - resp.match);
    
    resps = [resps; resp];
end

%% test plot
figure
hold on

% for legends
lbl1 = [];
lbl2 = [];
lbl3 = [];

% loop over k's
s = size(stim.ks,2);
for i=1:s
    % match
    resps(i).match_corr_r1 = mean(resps(i).match_corr(:,end-1:end), 2);
    resps(i).match_corr_r0 = mean(resps(i).match_corr(:,1:end-2), 2);
    resps(i).match_corr_mean = mean([resps(i).match_corr_r1, resps(i).match_corr_r0],2);
    
    subplot(2,2,1);hold on;
    plot(stim.eps, resps(i).match_corr, '--', 'Color',[i/s, 0, 1-(i/s)], 'LineWidth', 1);
    lbl1 = [lbl1 plot(stim.eps, resps(i).match_corr_mean, 'Color',[i/s, 0, 1-(i/s)], 'LineWidth', 4)];

    % center
    resps(i).center_corr_r1 = mean(resps(i).center_corr(:,end-1:end), 2);
    resps(i).center_corr_r0 = mean(resps(i).center_corr(:,1:end-2), 2);
    resps(i).center_corr_mean = mean([resps(i).center_corr_r1, resps(i).center_corr_r0],2);
    
    subplot(2,2,3);hold on;
    plot(stim.eps, resps(i).center_corr, '--', 'Color',[i/s, 0.7, 1-(i/s)], 'LineWidth', 1);
    lbl2 = [lbl2 plot(stim.eps, resps(i).center_corr_mean, 'Color',[i/s,0.7, 1-(i/s)], 'LineWidth', 4)];
    
    % difference    
    resps(i).diff_corr = resps(i).match_corr-resps(i).center_corr;
    resps(i).diff_corr_mean = resps(i).match_corr_mean - resps(i).center_corr_mean;
    
    subplot(2,2,[2 4]);hold on;
    plot(stim.eps, resps(i).diff_corr, '--', 'Color',[i/s, 0.35, 1-(i/s)], 'LineWidth', 1);
    lbl3 = [lbl3 plot(stim.eps, resps(i).diff_corr_mean, 'Color',[i/s,0.35, 1-(i/s)], 'LineWidth', 4)];
end


subplot(2,2,1);
xlim(stim.eps_range);
title(sprintf("Psychometric: Match"));
xlabel("stim: dist from center");ylabel("P(correct)");
legend(lbl1,stim.kslab,'Location','best');

subplot(2,2,3);
xlim(stim.eps_range);
title(sprintf("Psychometric: Center"));
xlabel("stim: dist from center");ylabel("P(correct)");
legend(lbl2,stim.kslab,'Location','best');

subplot(2,2,[2 4]);
xlim(stim.eps_range);
title(sprintf("Difference in curves (M-C)"));
xlabel("stim: dist from center");ylabel("\Delta P(correct)");
legend(lbl3,stim.kslab,'Location','best');


%% max diff

for i=1:s
    diff = resps(i).match_corr-resps(i).center_corr;
    diff_mean = resps(i).match_corr_mean - resps(i).center_corr_mean;
    
    f = figure;
    clim = [0 1];
    
    ax1 = subplot(1,2,1);
    hold on
    
    title("match");
    
    imagesc(resps(i).match_corr, clim);
    ylabel("P(correct)");
    xlabel("W stimuli");
    
    ax2 = subplot(1,2,2);
    hold on
    
    imagesc(resps(i).center_corr, clim);
    ylabel("P(correct)")
    xlabel("W stimuli")
    colorbar;
end







