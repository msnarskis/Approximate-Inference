%% Preprocess and Parameters
clear all;
clc;

% Stimulus Parameters
eps_range = [-4, 4];
num_eps = 100;
eps = linspace(eps_range(1), eps_range(2), num_eps);

k = 2;
w = 1;

dmu = (range(2) - range(1)) / range(3);

% Model Parameters
pr_R = 0.5;
% num_pr_R = 3;
% pr_R = linspace(0, 1, num_pr_R);
num_pr_C = 1;
pr_C = linspace(0, 1, num_pr_C);
num_pr_W = 1;
pr_W = linspace(0.8, 1, num_pr_W);

sig_v = 0.05;
sig_t = 0.5;
sig_n = 0.5;
sig_ap = 100;
sig_vp = 100;
mu_vp = 5;

lapse = 0;
nsamples = 1;

%% Generate Stimulus (k, eps, ext_smp)

ext_samp = 10000; % draws of X ~ eps

mu_hat = eps / ((2*w - 1));

X_t_nW = randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_t) + repmat(mu_hat, [k, 1, ext_samp]);
X_n_nW = randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_n) - repmat(mu_hat, [k, 1, ext_samp]);
X_vr1 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v) + abs(repmat(mu_hat, [k, 1, ext_samp])); % audio cue locations for each k (matched)
X_vl1 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v) - abs(repmat(mu_hat, [k, 1, ext_samp])); % audio cue locations for each k (matched)
X_vr0 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v); % audio cue locations for each k (central)
X_vl0 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v); % audio cue locations for each k (central)

W = sign(binornd(1, w, size(X_t_nW)) - 0.5); % correct response
X_t = W .* X_t_nW;
X_n = W .* X_n_nW;

% TEST: check mean(X_t, 3) appraoches eps [YES]

%% Model Predictions

results_c1 = zeros([num_eps, num_pr_C, num_pr_W]); % matched
results_c0 = zeros([num_eps, num_pr_C, num_pr_W]); % central

for i = 1:num_pr_C
    for j = 1:num_pr_W
        
        theta = [pr_R, pr_C(i), pr_W(j), sig_v, sig_t, sig_n, sig_ap, sig_vp, mu_vp, lapse, nsamples];
        
        results_c1(:, i, j) = generate_psychometric_curve_2(theta, X_t, X_n, X_vr1, X_vl1); % matched
        results_c0(:, i, j) = generate_psychometric_curve_2(theta, X_t, X_n, X_vr0, X_vl0); % central
        
    end % pr_W
end % pr_C


%% Analysis & Plot
figure(1)
hold on

for i = 1:num_pr_C
    for j = 1:num_pr_W
        
        plot(mu_hat, results_c1(:,num_pr_C,j), ':', 'linewidth', 0.2, 'color', 'r');
        plot(mu_hat, results_c0(:,i,j), ':', 'linewidth', 0.2, 'color', 'b');
        
    end % pr_W
end % pr_C

c1 = mean(mean(results_c1 ,3), 2); % matched
c0 = mean(mean(results_c0 ,3), 2); % central

mtch = plot(mu_hat, c1, 'linewidth', 2, 'color', 'r');
cnt = plot(mu_hat, c0, 'linewidth', 2, 'color', 'b');

%% diff
c1 = mean(mean(results_c1 ,3), 2); % matched
c0 = mean(mean(results_c0 ,3), 2); % central

%figure(2)
hold on

%diff = plot(mu_hat, (c1-c0)+0.5, 'linewidth', 2);

legend([mtch, cnt], {"matched", "central"}, 'Location', 'best');
