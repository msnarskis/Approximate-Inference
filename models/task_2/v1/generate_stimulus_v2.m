%generate_stimulus_v2([0, 5], 10, 3, [1 1 1; 1 -1 1; -1 -1 1; -1 -1 -1])
function [stim] = generate_stimulus_v2(eps_range, n, k, W)
    % Generates (not noisy) Stimulus (e_t) for 
    % corresponding values of W over n stimulus
    % size(stim) = (k: values, n: location, W: corespond to R)

    loc = linspace(eps_range(1), eps_range(2), n); % abs val of locations
    loc_k = repmat(loc, k, 1); % repeated
    
    stim = zeros(k, n, size(W, 1)); % init
    
    for w = 1:size(W,1) % loop over W values
        for i = 1:n     % loop over k values
            stim(:,i,w) = loc_k(:,i) .* W(w,:)';
        end
    end
    
    % stim % DEBUG
end