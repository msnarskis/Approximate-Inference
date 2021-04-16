%% Stimulus Parameters
eps_range = [0, 2];
num_eps = 200;

k = 1;
w = 1;

dmu = (range(2) - range(1)) / range(3);

%% Model Parameters
pr_R = 0.5;
num_pr_C = 1;
pr_C = linspace(0, .5, num_pr_C);
num_pr_W = 1;
pr_W = linspace(0.5, .5, num_pr_W);

sig_v = 0.1;
sig_t = 1;
sig_n = 1;
sig_ap = 100;
sig_vp = 100;
mu_vp = 0;

lapse = 0;
nsamples = 1;

%% Model Pedictions
figure(4)
hold on

k = linspace(1, 10, 10);
w = linspace(0.5, 1, 10);
max_diff = zeros([10, 10, 2]);

for i=1:10
    for j=1:10  
        eps_range = [-2*sig_t/sqrt(k(i)), 2*sig_t/sqrt(k(i))];
        dmu = (eps_range(2) - eps_range(1))/num_eps;
                
        stim = [eps_range num_eps k(i) w(j)];
        model = [pr_R pr_C pr_W sig_v sig_t sig_n sig_ap sig_vp mu_vp lapse nsamples];

        [mu_hat, c0, c1] = approx_model(stim, model);
        
        plot(abs(c1-c0));
                
        max_diff(i, j, 1) = max(abs(c1-c0));
        max_diff(i, j, 2) = sum(abs(c1-c0))*dmu;
        max_diff(:,:,1)
    end
end

%% Plot

figure(3)
imagesc(w, k, max_diff(:,:,1));
colorbar;
title("Max Difference in performance between Matched v Central");
xlabel("w (frequency stim on correct side)");
ylabel("k (number of frames)");

figure(2)
imagesc(w, k, max_diff(:,:,2));
colorbar;
title("plot1: sum abs diff m vs c");
xlabel("w");
ylabel("k");



