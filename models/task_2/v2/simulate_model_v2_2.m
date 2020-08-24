%% Model (v2.2) Testing

% stimulus parameters
eps_range = [0, 30];
n = 20;
k = 5;

% model parameters
sig_t = 10; % (std dev)
sig_n = 10;
sig_v = 1;
sig_sa = 100;
sig_sv = 100;

pr_R = 0.5;
pr_C = 0.5;

nsamp = 1; % trials per W vector
ntrials = 4000;

% noisy stimulus generation
noise_a = sig_t;
noise_v = sig_v;

%% Generate Stimulus
% for R = 0
W = ones(2^k,k);
for w = 1:2^k
    c = dec2bin(w-1,k);
    for i=1:k
        W(w,i) = W(w,i) * sign(str2num(c(i))-0.5);
    end
end

% correct answer
corr = repmat([ones(2^k-2,1);0;0],1,n)';

% TEST
W = [ones(1,k); zeros(1,k)];
corr = ones(n,2);

% stimulus var
[stim, W] = generate_stim_v2_2(eps_range, n, k);
eps = linspace(eps_range(1), eps_range(2), n);

%% Run Model

% returns resp: (n: locations, w: size(W))

center = model_v2_2(stim, 0, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp, ntrials);

%center_corr = abs(corr - center);

match = model_v2_2(stim, 1, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp, ntrials);

%match_corr = abs(corr - match);

%% Model Plotting
plot_debug_corr(eps, match, center,...
    sprintf('M(r) v C(b): k=%d, p(C)=%.1f nsamp=%d', k, pr_C, nsamp));

%% PLOT
figure
hold on
plot(eps, mean(match,2), 'r', 'LineWidth', 4);
plot(eps, mean(center,2), 'b', 'LineWidth', 4);
legend('matched','center')

