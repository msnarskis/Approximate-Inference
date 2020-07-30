% Simulate for R Bias
% Test to see how choice bias with no evidence
% changes as function of frames (k)

% stimulus parameters
K = 5;
ns = 20; % number of noise parameters
noise = linspace(0.1,10,ns);


% model parameters
sig_t = 10; % (std dev)
sig_n = 10;
sig_v = 1;
sig_sa = 100;
sig_sv = 100;

pr_R = 0.5;
pr_C = 0.5;

nsamp = 1; % trials per W vector
ntrials = 1000;

% noisy stimulus generation
noise_a = sig_t;
noise_v = sig_v;

%% Generate Stimulus and Simulate Models

% init vars
center = zeros(K-1,ns);
match = zeros(K-1,ns);

frames = 2:K;

for k=2:K
    [stim, w] = generate_stim_v2_2([0,0],1,k);
    k
    
    for n=1:ns
        
        %sig_t = noise(n);
        %sig_n = noise(n);
        %sig_v = noise(n)/10;
        noise_a = noise(n);
        noise_v = noise(n)/10;
        
        center(k-1,n) = mean(model_v2_2(stim, 0, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
            pr_R, pr_C, noise_a, noise_v, nsamp, ntrials),1);

        match(k-1,n) = mean(model_v2_2(stim, 1, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
            pr_R, pr_C, noise_a, noise_v, nsamp, ntrials),1);
    end
end

%% Plot

figure
hold on


col = [linspace(1,0,ns)' linspace(0,1,ns)' linspace(0,1,ns)'];

title('P(R=1|eps=0) by k - const param');
xlabel('k');ylabel('P(R=1)');
for n=1:ns
    plot(frames, center(:,n), 'LineWidth', 3, 'color', col(n,:));
end
for n=1:ns
    plot(frames, match(:,n), '--', 'LineWidth', 3, 'color', col(n,:));
end




