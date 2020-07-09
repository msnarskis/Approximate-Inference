function y = generate_psychometric_curve(thet, epsilon_aud, epsilon_vis_left)

% thet - Parameters (passed as a row)
%   1 - prior probibility of right response,
%   2 - prior probability of same cause,
%   3 - visual uncertainty ($\sigma_v^2$)
%   4 - auditory uncertainty ($\sigma_a^2$)
%   5 - lapse rate
%   6 - nsamples (1 - 100), -Inf, or Inf

% epsilon_aud - Auditory Cue Location
% epsilon_vis_left - Left Visual Cue Location

% DONT FORGET: convert to response using binornd(y)!!!!!

ext_samp = 10000;

epsilon_aud = epsilon_aud(:)';
epsilon_vis_left = epsilon_vis_left(:)';

k = size(epsilon_aud, 2); % number of stimuli
n = size(thet, 1); % num parameters

epsilon_aud = repmat(epsilon_aud, [n, 1, ext_samp]);
epsilon_vis_left = repmat(epsilon_vis_left, [n, 1, ext_samp]);


pr1 = repmat(thet(:,1), [1, k, ext_samp]);
pc1 = repmat(thet(:,2), [1, k, ext_samp]);

sig_v_eff_2 = repmat(thet(:,3), [1 ,k, ext_samp]);
sig_a_eff_2 = repmat(thet(:,4), [1, k, ext_samp]);

lam = repmat(thet(:,5), [1, k, ext_samp]);
nsamp = round(repmat(thet(:,6), [1, k, ext_samp]));

sig_av_scl_2 = sig_a_eff_2 + sig_v_eff_2;
alp_av = sig_v_eff_2 ./ sig_av_scl_2;

x_a_eff = epsilon_aud + normrnd( zeros(size(epsilon_aud)), sqrt(sig_a_eff_2));
x_v_eff = epsilon_vis_left + normrnd( zeros(size(epsilon_aud)), sqrt(sig_v_eff_2));

sig_av_comb_2 = sig_a_eff_2 .* alp_av;


num = zeros([n, k, ext_samp]);
den = zeros([n, k, ext_samp]);
% pbseps_a = repmat(eps_a, [1, 1, ext_samp])

tmpvar = log(0.5 * (1 + erf(-x_v_eff./sqrt(2*sig_v_eff_2))));

num(:,:,:,1) = log(pr1) + log(1-pc1) ...
    + log(0.5 * (1 + erf(x_a_eff./sqrt(2*sig_a_eff_2)))) + tmpvar;

num(:,:,:,2) = log(pr1) + log(pc1) ...
    + log(normpdf(x_a_eff,-x_v_eff,sqrt((sig_av_scl_2)))) ...
    + log(0.5 + 0.5 * erf( (x_a_eff.*alp_av-x_v_eff.*(1-alp_av)) ./ sqrt(2*sig_av_comb_2)));

den(:,:,:,1) = num(:,:,:,1);
den(:,:,:,2) = num(:,:,:,2);
den(:,:,:,3) = log(1-pr1) + log(1-pc1) ...
    + log(0.5 * (1 + erf(-x_a_eff./sqrt(2*sig_a_eff_2)))) + tmpvar;
den(:,:,:,4) = log(1-pr1) + log(pc1) ...
    + log(normpdf(x_a_eff, x_v_eff, sqrt((sig_av_scl_2)))) ...
    + log(0.5 + 0.5 * erf( -1*(x_a_eff.*alp_av + x_v_eff.*(1-alp_av)) ./ sqrt(2*sig_av_comb_2)));

lp=logsumexp(num,4)-logsumexp(den,4);
prb_right=exp(lp);

prb_right(isnan(prb_right))=0.5;

pb=zeros(size(prb_right));
pb(nsamp~=Inf)=binocdf(floor(nsamp(nsamp~=Inf)/2),nsamp(nsamp~=Inf),prb_right(nsamp~=Inf),'upper')+0.5*binopdf(nsamp(nsamp~=Inf)/2,nsamp(nsamp~=Inf),prb_right(nsamp~=Inf));
pb(nsamp==Inf)=(prb_right(nsamp==Inf)>0.5)+0.5*(prb_right(nsamp==Inf)==0.5);
%pb = prb_right;
phat=0.5.*lam+(1-lam).*(pb);

y=mean(phat,3);

end