%% Stimulus Parameters
eps_range = [0, 5];
num_eps = 50;

num_k = 3;
K = 1:2:(2*num_k-1);
num_w = 5;
W = linspace(0.5, 1, num_w);

resp0 = zeros([num_k, num_w, num_eps]);
resp1 = zeros([num_k, num_w, num_eps]);
diff = zeros([num_k, num_w, num_eps]);
time = [];

%% Model Parameters
pr_R = 0.5;
pr_C = 0.5;
pr_W = 1;

sig_v = 1;
sig_t = 10;
sig_n = 10;
sig_ap = 100;
sig_vp = 100;
mu_vp = 0;

lapse = 0;
nsamples = 1;

%% Model Pedictions
figure;
hold on;

for k=K
    for w=W
		sprintf("k=%d\tw=%1.2f ",k,w) % UX GUI
		sTime = cputime;

		% set params
        stim = [eps_range, num_eps, k, w];
        model = [pr_R pr_C 1 sig_v sig_t sig_n sig_ap sig_vp mu_vp lapse nsamples];

        [mu_hat, c0, c1] = approx_model(stim, model);
		
		time = [time, (cputime - sTime)] % sim done timestamp
		
        plot(mu_hat, c0, 'LineWidth', 2, 'DisplayName', sprintf("c:   k=%d ",k,w));
		plot(mu_hat, c1, 'LineWidth', 2, 'DisplayName', sprintf("m: k=%d ",k,w));
    end
end

legend('Location','eastoutside');

%% Plot
% 
% figure
% imagesc(w, k, max_diff(:,:,1));
% colorbar;
% title("Max Difference in performance between Matched v Central");
% xlabel("w (frequency stim on correct side)");
% ylabel("k (number of frames)");
% 
% figure
% imagesc(w, k, max_diff(:,:,2));
% colorbar;
% title("plot1: sum abs diff m vs c");
% xlabel("w");
% ylabel("k");



