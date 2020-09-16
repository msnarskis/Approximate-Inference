%% Model (v2.2) Testing
clc, clear all;

% stimulus and plot selection
gen_stim = @(eps_range, n, k) generate_stim_v2_2(eps_range, n, k);
% gen_stim = @(eps_range, n, k) gen_stim_test_bias(k);

plot_me = @(eps, match_corr, center_corr, k, pr_C, nsamp, W) ...
    plot_debug_corr(eps, match_corr, center_corr, ...
    sprintf('M(r) v C(b): k=%d, p(C)=%.1f nsamp=%d', k, pr_C, nsamp));
% plot_me = @(eps, match_corr, center_corr, k, pr_C, nsamp, W) ...
%     plot_test_bias(match_corr, center_corr, W);

% stimulus parameters
eps_range = [0, 30];
n = 20;
k = 3;

% model parameters
sig_t = 10; % (std dev)
sig_n = 10;
sig_v = 1;
sig_sa = 100;
sig_sv = 100;

pr_R = 0.5;
pr_C = 0.5;

nsamp = 1; % trials per W vector
ntrials = 2000;

% noisy stimulus generation
noise_a = sqrt(sig_t);
noise_v = sqrt(sig_v);


%% Generate Stimulus

% stimulus var
[stim, W, corr] = gen_stim(eps_range, n, k);
eps = linspace(eps_range(1), eps_range(2), n);

stim = stim;


%% Run Model

% returns resp: (n: locations, w: size(W))

center = model_v2_2(stim, 0, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp, ntrials);

match = model_v2_2(stim, 1, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp, ntrials);

match_corr = abs(corr - match);
center_corr = abs(corr - center);


%% Model Plotting

plot_me(eps, match_corr, center_corr, k, pr_C, nsamp, W);

