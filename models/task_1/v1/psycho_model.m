%% Generate psychometric curve for given parameter values

function [mu_a, c1, c0] = psycho_model(theta, stim)

%Use theta and stimulus parameters to generate response curve
% 
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
%   11 - nsamples (1 - 100), -Inf, or Inf
%
% stim - Parameters
%   1, 2 - range: audio cue loc beginning and end
%   3 - number of locations to look at
%   4 - sig_t: variance of audio location around mean
%   5 - sig_n
%   6 - sig_v
%   7 - k number of stimuli shown in succesion
%   8 - probability that audio cue is at mu_a_hat vs -mu_a_hat
%       where mu_a_hat is adjusted to p s.t. mean audio cue location
%       is at mu_a
%   9 - number of samples to average over
%
% eps_a  - Auditory Cue Location     size(num_locations, k)
% eps_vr - Left Visual Cue Location size(num_locations, k)
%
% match - 1: matched condition (eps_v = eps_a)
%         0: center condition  (eps_v = zeros)

% generate all X_*i independantly

mu_a = linspace(stim(1), stim(2), stim(3));

[eps_t, eps_n, eps_vr1, eps_vl1, eps_vr0, eps_vl0] = generate_stimulus(mu_a, stim(4), stim(5), stim(6), stim(7), stim(8), stim(9));

c1 = generate_psychometric_curve_2(theta, eps_t, eps_n, eps_vr1, eps_vl1);
c0 = generate_psychometric_curve_2(theta, eps_t, eps_n, eps_vr0, eps_vl0);

end