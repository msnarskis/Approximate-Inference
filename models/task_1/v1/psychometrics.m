%% Generate psychometric curves for various conditions

clear all;
clc;

% stim - Parameters
%   1, 2 - range: audio cue loc beginning and end
%   3 - number of locations to look at
%   4 - sig_a: variance of audio location around mean
%   5 - k number of stimuli shown in succesion
%   6 - probability that audio cue is at mu_a_hat vs -mu_a_hat
%       where mu_a_hat is adjusted to p s.t. mean audio cue location
%       is at mu_a
%   7 - number of samples to average over (for stim_sig_a > 0)
%
% eps_a  - Auditory Cue Location     size(num_locations, k)
% eps_vr - Left Visual Cue Location  size(num_locations, k)
%
% match - 1: matched condition (eps_v = eps_a)
%         0: center condition  (eps_v = zeros)

% stimulus variables
range = [0, 3, 600];
stim_sig_t = 1;
stim_sig_n = 1;
stim_sig_v = .1;
noise = [stim_sig_t, stim_sig_n, stim_sig_v];
k = 1;
w = 1;
ext_samp = 10000;

stim = [range, noise, k, w, ext_samp];

dmu = (range(2) - range(1)) / range(3);

% thet - Parameters (passed as a row)
%   1  - prior probibility of right response  P(R)
%   2  - prior probability of same cause      P(C)
%   3  - prior probability of wrong side      P(W)
%   4  - visual uncertainty                   \sigma_v^2    
%   5  - tone uncertainty                     \sigma_t^2
%   6  - noise uncertainty                    \sigma_n^2
%   7  - auditory perception uncertainty      \sigma_ap^2
%   8  - visual perception uncertainty        \sigma_vp^2
%   9  - visual perception prior mean         \mu_vp^2
%   10 - lapse rate
%   11 - nsamples (1 - 100), -Inf, or Inf5

% theta variables
pr_R = .5; 
pr_C = .5;
pr_W = 1;
sig_v = stim_sig_v;
sig_t = stim_sig_t;
sig_n = stim_sig_n;
sig_ap = 100;
sig_vp = 100;
mu_vp = 0;
lapse = 0;
nsamples = 1;

theta = [pr_R, pr_C, pr_W, sig_v, sig_t, sig_n, sig_ap, sig_vp, mu_vp, lapse, nsamples];


% generate values
theta(11) = 1; % nsamp
[mu_a1, c1_1, c0_1] = psycho_model(theta, stim);

%theta(11) = 100; % nsamp
%[mu_a2, c1_2, c0_2] = psycho_model(theta, stim);

%theat(11) = Inf; % k
%[mu_a3, c1_3, c0_3] = psycho_model(theta, stim);

stim(5) = 4; % k
%[mu_a4, c1_4, c0_4] = psycho_model(theta, stim);

stim(5) = 5; % k
%[mu_a5, c1_5, c0_5] = psycho_model(theta, stim);


% diff between match and central
diff_1 = c1_1 - c0_1;
%diff_2 = sum(abs(c1_2 - c0_2)) * dmu;
%diff_3 = sum(abs(c1_3 - c0_3)) * dmu;
%diff_4 = sum(abs(c1_4 - c0_4)) * dmu;
%diff_5 = sum(abs(c1_5 - c0_5)) * dmu;

%% plot(50:end)
figure;
hold on
%plot([1 2 3 4 5], [diff_1 diff_2 diff_3 diff_4 diff_5])
plot(mu_a1, c1_1, 'Linewidth', 2);
plot(mu_a1, c0_1, 'Linewidth', 2);
%plot(mu_a2, c1_2, 'Linewidth', 2);
%plot(mu_a2, c0_2, 'Linewidth', 2);
%plot(mu_a3, c1_3, 'Linewidth', 2);
%plot(mu_a3, c0_3, 'Linewidth', 2);
%plot(mu_a5, c0_5, 'Linewidth', 2);
%plot(mu_a1, .5*ones(size(mu_a1)), 'k--');
plot([0,0], [0,1], 'k--');
legend('match smp=1', 'central smp=1', 'match smp=100', 'central smp=100', 'match smp=\infty', 'central smp=\infty', 'location', 'best')
title("P(R|X) as a function of k");

%% TEST
figure;
hold on
plot(movmean(c0_1,12), movmean(diff_1,12))

