%% Model (v2.2) Testing
clc, clear all;

% stimulus and plot selection
% gen_stim = @(stim) generate_stim_v2_2(stim);
gen_stim = @(stim) gen_stim_test_bias(stim);

% plot_me = @(stim,resp,par) ...
%     plot_debug_corr(stim,resp, ...
%     sprintf('M(r) v C(b): k=%d, p(C)=%.1f nsamp=%d', stim.k, par.pr_C, par.nsamp));
plot_me = @(stim,resp,par) ...
    plot_test_bias(stim,resp);

% stimulus parameters
stim.eps_range = [0, 10];
stim.n = 20;
stim.k = 3;

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


%% Generate Stimulus

stim = gen_stim(stim);

%% Run Model

% returns resp: (n: locations, w: size(W))

resp.center = model_v2_2(stim, par, 0);

resp.match = model_v2_2(stim, par, 1);

resp.match_corr = abs(stim.corr - resp.match);
resp.center_corr = abs(stim.corr - resp.center);


%% Model Plotting

plot_me(stim, resp, par);

