function [mu_hat results_c0 results_c1] = approx_model(stim, model)

    % Stimulus Parameters
    eps_range = [stim(1) stim(2)];
    num_eps = stim(3);
    eps = linspace(eps_range(1), eps_range(2), num_eps);

    k = stim(4);
    w = stim(5);

    dmu = (range(2) - range(1)) / range(3);

    % Model Parameters
    sig_v = model(4);
    sig_t = model(5);
    sig_n = model(6);

    %% Generate Stimulus
    
    ext_samp = 5000; % draws of X ~ eps

    mu_hat = eps;% / ((2*w - 1));

    X_t_nW = randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_t) + repmat(mu_hat, [k, 1, ext_samp]);
    X_n_nW = randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_n) - repmat(mu_hat, [k, 1, ext_samp]);
    X_vr1 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v) + abs(repmat(mu_hat, [k, 1, ext_samp])); % audio cue locations for each k (matched)
    X_vl1 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v) - abs(repmat(mu_hat, [k, 1, ext_samp])); % audio cue locations for each k (matched)
    X_vr0 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v); % audio cue locations for each k (central)
    X_vl0 =  randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_v); % audio cue locations for each k (central)

    W = sign(binornd(1, w, size(X_t_nW)) - 0.5); % correct response
    X_t = W .* X_t_nW;
    X_n = W .* X_n_nW;
    
    %% Model Predictions

    results_c1 = generate_psychometric_curve_2(model, X_t, X_n, X_vr1, X_vl1); % matched
    results_c0 = generate_psychometric_curve_2(model, X_t, X_n, X_vr0, X_vl0); % central


end