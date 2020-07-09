function [eps_t, eps_n, eps_vr1, eps_vl1, eps_vr0, eps_vl0] = generate_stimulus(mu_a, sig_t, sig_n, sig_v, k, w, ext_samp)

%Generate stimulus: row vector of length k where eps_a = eps_v
% and ~ w * N(mu_ahat, sig_a) + (1-w) * N(-mu_ahat, sig_a)
% where mu_ahat is s.t. the mean of eps_a is mu_a
%
% condition 0  -->  p = 1,       k = 1
% condition 1  -->  p = 1,       k < 1
% condition 2  -->  p in (.5,1], k < 1
%
% noise:    1 - sig_t
%           2 - sig_nmu_a / ((2*w - 1)
%           3 - sig_v

mu_hat = mu_a / ((eps + 2*w - 1)); % *sqrt(k)); % effective mean

eps_t_nW = randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_t) + repmat(mu_hat, [k, 1, ext_samp]); 
eps_n_nW = randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_n) - repmat(mu_hat, [k, 1, ext_samp]); 
eps_vr1 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v) + abs(repmat(mu_hat, [k, 1, ext_samp])); % audio cue locations for each k
eps_vl1 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v) - abs(repmat(mu_hat, [k, 1, ext_samp])); % audio cue locations for each k
eps_vr0 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v); % audio cue locations for each k
eps_vl0 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v); % audio cue locations for each k

W = sign(binornd(1, w, size(eps_t_nW)) - 0.5); % correct response
sum(count(string(W), '-1'), 2);

eps_t = W .* eps_t_nW;
eps_n = W .* eps_n_nW;
end





