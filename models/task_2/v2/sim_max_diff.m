%% Simulating for max_diff between conditions
clc; clear all;

% stimulus parameters
stim.eps_range = [0, 8];
stim.n = 30;
stim.ks = [2,3,5,7];
stim.kslab = {'k=2','k=3','k=5','k=7'};

s = size(stim.ks,2);

% model parameters
par.var_t = 10; % (std dev)
par.var_n = 10;
par.var_v = 1;
par.var_sa = 100;
par.var_sv = 100;

par.pr_R = 0.5;
par.pr_C = 0.5;

par.nsamp = 100; % trials per W vector
par.ntrials = 1000;

par.noisy_in = 1;

%% Simulate
model = Model_Task_2(par);

resps = [];

time = [];

for k=stim.ks
    sprintf("k=%d",k) % for keeping track during run
    sTime = cputime;
    
    stim.k = k;
    resp.k = k;
    
    stim = generate_stim_v2_2(stim);
    
    resp.center = model.simulate(stim,0);
    resp.match  = model.simulate(stim,1);
    
    resp.center_corr = abs(stim.corr - resp.center);
    resp.match_corr  = abs(stim.corr - resp.match);
    
    resps = [resps; resp];
    
    time = [time, (cputime - sTime)]
end

%% plot psychometric curves and diff between match, center
figure
hold on

% for legends
lbl1 = [];
lbl2 = [];
lbl3 = [];

% loop over k's
for i=1:s
    % match
    resps(i).match_corr_r1 = mean(resps(i).match_corr(:,end-1:end), 2);
    resps(i).match_corr_r0 = mean(resps(i).match_corr(:,1:end-2), 2);
    resps(i).match_corr_mean = mean([resps(i).match_corr_r1, resps(i).match_corr_r0],2);
    
    subplot(2,2,1);hold on;
    plot(stim.eps, resps(i).match_corr, '--', 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 1);
    lbl1 = [lbl1 plot(stim.eps, resps(i).match_corr_mean, 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 5)];

    % center
    resps(i).center_corr_r1 = mean(resps(i).center_corr(:,end-1:end), 2);
    resps(i).center_corr_r0 = mean(resps(i).center_corr(:,1:end-2), 2);
    resps(i).center_corr_mean = mean([resps(i).center_corr_r1, resps(i).center_corr_r0],2);
    
    subplot(2,2,3);hold on;
    plot(stim.eps, resps(i).center_corr, '--', 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 1);
    lbl2 = [lbl2 plot(stim.eps, resps(i).center_corr_mean, 'Color',[i/s,0.3, 1-(i/s)], 'LineWidth', 5)];
    
    % difference    
    resps(i).diff_corr = resps(i).match_corr-resps(i).center_corr;
    resps(i).diff_corr_mean = resps(i).match_corr_mean - resps(i).center_corr_mean;
    
    subplot(2,2,[2 4]);hold on;
    plot(stim.eps, resps(i).diff_corr, '--', 'Color',[i/s, 0.3, 1-(i/s)], 'LineWidth', 1);
    lbl3 = [lbl3 plot(stim.eps, resps(i).diff_corr_mean, 'Color',[i/s,0.3, 1-(i/s)], 'LineWidth', 5)];
end


subplot(2,2,1);
xlim(stim.eps_range);
title(sprintf("Psychometric: Match -- nsamp:%d",par.nsamp));
xlabel("stim: dist from center");ylabel("P(correct)");
legend(lbl1,stim.kslab,'Location','best');

subplot(2,2,3);
xlim(stim.eps_range);
title(sprintf("Psychometric: Center -- nsamp:%d",par.nsamp));
xlabel("stim: dist from center");ylabel("P(correct)");
legend(lbl2,stim.kslab,'Location','best');

subplot(2,2,[2 4]);
xlim(stim.eps_range);
title(sprintf("Difference in curves (M-C) -- nsamp:%d",par.nsamp));
xlabel("stim: dist from center");ylabel("\Delta P(correct)");
legend(lbl3,stim.kslab,'Location','best');

set(gcf,'Position',[100 100 1600 1000]);


%% plot flattened curves max diff

for i=1:s
    f = figure;
    clim = [0 1];
    xlims = [0.5, size(resps(i).match_corr,1)+0.5];
    ylims = [0.5, size(resps(i).match_corr,2)+0.5];
    
    % match condition
    ax1 = subplot(2,2,1);
    hold on
    
    imagesc(resps(i).match_corr', clim);
    xlabel("stimulus location");ylabel("stimulus class");
    xlim(xlims);ylim(ylims);yticks([]);
    title(sprintf("Flattened psychometric curves: Match; k=%d, nsamp=%d",stim.ks(i),par.nsamp));
    c = colorbar; ylabel(c, "P(correct)");
    
    % center condition
    ax2 = subplot(2,2,3);
    hold on

    imagesc(resps(i).center_corr', clim);
    xlabel("stimulus location");ylabel("stimulus class");
    xlim(xlims);ylim(ylims);yticks([]);
    title(sprintf("Flattened psychometric curves: Center; k=%d, nsamp=%d",stim.ks(i),par.nsamp));
    c = colorbar; ylabel(c, "P(correct)");
        
    % difference between condsd
    ax3 = subplot(2,2,[2,4]);
    hold on
    
    imagesc(resps(i).diff_corr');
    xlabel("stimulus location");ylabel("stimulus class");
    xlim(xlims);ylim(ylims);yticks([]);
    title(sprintf("Difference btw Match - Center; k=%d, nsamp=%d",stim.ks(i),par.nsamp));
    c = colorbar; ylabel(c, "\Delta P(correct)");
    
    set(gcf,'Position',[100 100 1600 1000]);
end

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

%% Save
file = strcat(pwd, sprintf("/models/task_2/v2/sims/sim_max_diff-t:%d-k:%d-nsamp:%d.mat", par.ntrials, stim.ks(end), par.nsamp));
save(file,'resps','par','stim');





