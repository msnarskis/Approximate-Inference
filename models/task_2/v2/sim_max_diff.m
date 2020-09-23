%% Simulating for max_diff between conditions
clc, clear all;

% stimulus parameters
stim.eps_range = [0, 8];
stim.n = 30;
stim.ks = [2,3,5,7];

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

s = size(stim.ks,2);
for i=1:s
    % resp = mean([ mean(resp(i).match_corr(:,end-1:end),2), mean(resp.match_corr(:,1:end-2),2) ],2);
    resps(i).match_corr_r1 = mean(resps(i).match_corr(:,end-1:end), 2);
    resps(i).match_corr_r0 = mean(resps(i).match_corr(:,1:end-2), 2);
    resps(i).match_corr_mean = mean([resps(i).match_corr_r1, resps(i).match_corr_r0],2);
    
    plot(stim.eps, resps(i).match_corr, ':', 'Color',[i/s, 0, 1-(i/s)]);
    plot(stim.eps, resps(i).match_corr_mean, 'Color',[i/s, 0, 1-(i/s)], 'LineWidth', 3);

    resps(i).center_corr_r1 = mean(resps(i).center_corr(:,end-1:end), 2);
    resps(i).center_corr_r0 = mean(resps(i).center_corr(:,1:end-2), 2);
    resps(i).center_corr_mean = mean([resps(i).center_corr_r1, resps(i).center_corr_r0],2);
    
    plot(stim.eps, resps(i).center_corr, ':', 'Color',[i/s, 0.7, 1-(i/s)]);
    plot(stim.eps, resps(i).center_corr_mean, 'Color',[i/s,0.7, 1-(i/s)], 'LineWidth', 3);
end

%% max diff

for i=1:s
    diff = resps(i).match_corr-resps(i).center_corr;
    
    figure;
    title(sprintf("max diff m-c K=%d",stim.ks(i)));
    imagesc(diff);
    ylabel("max \Delta P(correct)")
    xlabel("W stimuli")
    colorbar;
end