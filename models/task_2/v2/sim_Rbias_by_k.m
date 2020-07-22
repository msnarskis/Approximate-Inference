% Simulate for R Bias
% Test to see how choice bias with no evidence
% changes as function of frames (k)

% stimulus parameters
K = 6;

% model parameters
sig_t = 10; % (std dev)
sig_n = 10;
sig_v = 1;
sig_sa = 100;
sig_sv = 100;

pr_R = 0.5;
pr_C = 0.5;

nsamp = 1; % trials per W vector
ntrials = 3000;

% noisy stimulus generation
noise_a = sig_t;
noise_v = sig_v;

%% Generate Stimulus and Simulate Models

% init vars
center = zeros(K-1,1);
match = zeros(K-1,1);

frames = 2:K;

for k=2:K
    [stim, w] = generate_stim_v2_2([0,0],1,k);
 
    center(k-1) = mean(model_v2_2(stim, 0, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
        pr_R, pr_C, noise_a, noise_v, nsamp, ntrials),1);

    match(k-1) = mean(model_v2_2(stim, 1, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
        pr_R, pr_C, noise_a, noise_v, nsamp, ntrials),1);
    
end

%% Plot

figure
hold on

title('P(R=1|eps=0) by k - 1');
xlabel('k');ylabel('P(R=1)');
plot(frames, center, 'r', 'LineWidth', 2);
plot(frames, match, 'b', 'LineWidth', 2);
legend('center','matched');
