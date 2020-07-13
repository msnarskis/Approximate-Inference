%% Model (v2.2) Testing

% stimulus parameters
eps_range = [0, 30];
k = 2;
n = 4;

% model parameters
sig_t = 10; % (std dev)
sig_n = 10;
sig_v = 1;
sig_sa = 100;
sig_sv = 100;

pr_R = 0.5;
pr_C = 0.5;

nsamp = 1; % trials per W vector
ntrials = 1;

% noisy stimulus generation
noise_a = 0;
noise_v = 0;

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
corr = repmat([0;ones(2^k-2,1);0],1,n)';

% stimulus var
stim = generate_stimulus_v2(eps_range, n, k, W);
eps = linspace(eps_range(1), eps_range(2), n);

%% Run Model

% returns resp: (n: locations, w: size(W))

center = model_v2_2(stim, 0, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp, ntrials);

center_corr = abs(corr - center);

match = model_v2_2(stim, 1, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp, ntrials);

match_corr = abs(corr - match);

%% Model Analyses

% match per W
figure;
hold on;

title("Match R=1");
xlabel("stimulus location (eps)");
ylabel("prob R=1");
ylim([0 1]);
w=size(W,1);
for i = 2:(w-1)
    plot(eps, match_corr(:,i), 'LineWidth', 4);
end

for i = [1,w]
    plot(eps, match_corr(:,i), '--', 'LineWidth', 4);
end


% center per W
figure;
hold on;

title("Center R=1");
xlabel("stimulus location (eps)");
ylabel("prob R=1");
ylim([0 1]);

for i = 2:(w-1)
    plot(eps, center_corr(:,i), 'LineWidth', 4);
end

for i = [1,w]
    plot(eps, center_corr(:,i), '--', 'LineWidth', 4);
end

%% 2D of above
% figure
% hold on
%
% title("Match")
% imagesc(eps, w_amp, match');
% xlabel("eps")
% ylabel("W")
% colorbar;
%
% figure
% hold on
%
% title("Center")
% imagesc(eps, w_amp, center');
% xlabel("eps")
% ylabel("W")
% colorbar;

%% Match vs Center
figure
hold on

title("Center (b) v Matched (r)");
xlabel("stimulus location (eps)");
ylabel("prob correct");
ylim([0 1]);

% for i = 1:w-2
%     plot(eps, match(:,i), 'r', 'LineWidth', 2);
%     plot(eps, center(:,i), 'b', 'LineWidth', 2);
% end
%
% for i = w-1:w
%     plot(eps, match(:,i), 'r--', 'LineWidth', 2);
%     plot(eps, center(:,i), 'b--', 'LineWidth', 2);
% end

plot(eps, mean(match_corr, 2), 'r', 'LineWidth', 4);
plot(eps, mean(center_corr, 2), 'b', 'LineWidth', 4);
