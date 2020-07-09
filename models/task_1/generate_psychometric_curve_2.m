function [psycho_curve] = generate_psychometric_curve_2(thet, eps_t, eps_n, eps_vr, eps_vl)

% see https://www.overleaf.com/read/cjstrytjpgqf for details

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

% eps_a - Auditory Cue Location     size(num_locations, k, smp)
% eps_vr - Left Visual Cue Location size(num_locations, k, smp)


num_loc = size(eps_t, 2); % number of a/v locations
k = size(eps_t, 1); % number of frames
ext_samp = size(eps_t, 3); % number of samples


% repeat matrices of parameters
pr_R  = thet(1);
pr_C = thet(2);
pr_W = thet(3);
sig_v = thet(4);
sig_t = thet(5);
sig_n = thet(6);
sig_ap = thet(7);
sig_vp = thet(8);
mu_vp = thet(9);

lr = thet(10);
nsamp = thet(11);

zero = zeros(size(eps_t));


% combined variables
a_tn  = sig_n ./ (sig_t + sig_n);
mu_tn = eps_t.*a_tn - eps_n.*(1-a_tn);
mu_v  = 0.5*(eps_vr - eps_vl); 

a_sa  = sig_ap ./ (sig_t.*a_tn + sig_ap);
mu_sa = mu_tn.*a_sa;
a_sv  = (2*sig_vp) ./ (sig_v + 2*sig_vp);
mu_sv = mu_v.*a_sv + mu_vp.*(1-a_sv);

a_aw  = (sig_v.*a_sv) ./ (2*sig_t.*a_tn.*a_sa + sig_v.*a_sv);
sig_aw = (sig_t.*a_tn.*a_sa) .* a_aw;


% calculations with R = 1
R = 1;

mu_aw  = mu_sa.*a_aw - R*mu_sv.*(1-a_aw);
mu_anw = mu_sa.*a_aw + R*mu_sv.*(1-a_aw);

nWnC = (1-pr_W).*(1-pr_C).*normcdf(zero, R*mu_sa, sqrt(sig_t.*a_tn.*a_sa))...
    .*normcdf(zero, -1*mu_sv, sqrt(0.5*sig_v.*a_sv));

WnC = (pr_W).*(1-pr_C).*normcdf(zero, -1*R*mu_sa, sqrt(sig_t.*a_tn.*a_sa))...
    .*normcdf(zero, -1*mu_sv, sqrt(0.5*sig_v.*a_sv));

nWC = (1-pr_W).*(pr_C).*normpdf(mu_sa, -1*R*mu_sv, sqrt(sig_t.*a_tn.*a_sa + .5*sig_v.*a_sv))...
    .*normcdf(zero, R*mu_aw, sqrt(sig_aw));

WC = (pr_W).*(pr_C).*normpdf(mu_sa, R*mu_sv, sqrt(sig_t.*a_tn.*a_sa + .5*sig_v.*a_sv))...
    .*normcdf(zero, -1*R*mu_anw, sqrt(sig_aw));


pX_r1 = exp(sum(log(nWnC + WnC + nWC + WC), 1)) * pr_R;


% calculations with R = -1
R = -1;

mu_aw  = mu_sa.*a_aw - R*mu_sv.*(1-a_aw);
mu_anw = mu_sa.*a_aw + R*mu_sv.*(1-a_aw);

nWnC = (1-pr_W).*(1-pr_C).*normcdf(zero, R*mu_sa, sqrt(sig_t.*a_tn.*a_sa))...
    .*normcdf(zero, -1*mu_sv, sqrt(0.5*sig_v.*a_sv));

WnC = (pr_W).*(1-pr_C).*normcdf(zero, -1*R*mu_sa, sqrt(sig_t.*a_tn.*a_sa))...
    .*normcdf(zero, -1*mu_sv, sqrt(0.5*sig_v.*a_sv));

nWC = (1-pr_W).*(pr_C).*normpdf(mu_sa, -1*R*mu_sv, sqrt(sig_t.*a_tn.*a_sa + .5*sig_v.*a_sv))...
    .*normcdf(zero, R*mu_aw, sqrt(sig_aw));

WC = (pr_W).*(pr_C).*normpdf(mu_sa, R*mu_sv, sqrt(sig_t.*a_tn.*a_sa + .5*sig_v.*a_sv))...
    .*normcdf(zero, -1*R*mu_anw, sqrt(sig_aw));


pX_r2 = exp(sum(log(nWnC + WnC + nWC + WC), 1)) * (1-pr_R);


% normalize, sample, and return
pr_Xr = pX_r1 ./ (pX_r1 + pX_r2);

pr_Xr(isnan(pr_Xr)) = .5;

% non_Inf samples
if nsamp~=Inf
    pr_X = binocdf(floor(nsamp/2), nsamp, pr_Xr, 'upper')...;
    + .5*binopdf(nsamp/2, nsamp, pr_Xr);
end

% Inf samples
if nsamp==Inf
    pr_X = (pr_Xr > 0.5) + 0.5*(pr_Xr == 0.5);
end

pr_X = mean(pr_X, 3);
psycho_curve = 0.5*lr + (1-lr)*pr_X;
end






