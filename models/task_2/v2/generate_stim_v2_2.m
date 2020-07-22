%generate_stimulus_v2([0, 5], 10, 3, [1 1 1; 1 -1 1; -1 -1 1; -1 -1 -1])
function [stim, W] = generate_stimulus_v2_2(eps_range, n, k)
    % Generates (not noisy) Stimulus (e_t) for 
    % all values of W over n locations
    % size(stim) = (k: frames, n: location, W: corresponding order)
    
    % generate W vec
    W = ones(2^k,k);
    for w = 1:2^k
        c = dec2bin(w-1,k);
        for i=1:k
            W(w,i) = W(w,i) * sign(str2num(c(i))-0.5);
        end
    end
    
    [~, i] = sort(abs(sum(W,2)));
    W = W(i,:);
    
    
    loc = linspace(eps_range(1), eps_range(2), n); % abs val of locations
    loc_k = repmat(loc, k, 1); % repeated
    
    stim = zeros(k, n, size(W, 1)); % init return var
    
    for w = 1:size(W,1) % loop over W values
        for i = 1:n     % loop over k values
            stim(:,i,w) = loc_k(:,i) .* W(w,:)';
        end
    end
    
end