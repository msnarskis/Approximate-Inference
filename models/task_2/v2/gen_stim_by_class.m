%generate_stimulus_v2([0, 5], 10, 3, [1 1 1; 1 -1 1; -1 -1 1; -1 -1 -1])
function [s] = gen_stim_by_class(s)
    % Generates (not noisy) Stimulus (e_t) for 
    % all values of W over n locations
    % size(stim) = (k: frames, n: location, W: corresponding order)
    
    % generate W vec
    s.W = ones(ceil(s.k/2),s.k);
    for w = 0:floor(s.k/2)
        s.W(w+1,:) = [ones(1,s.k-w),-1*ones(1,w)];
    end   
    
    loc = linspace(s.eps_range(1), s.eps_range(2), s.n); % abs val of locations
    loc_k = repmat(loc, s.k, 1); % repeated
    
    s.in = zeros(s.k, s.n, size(s.W, 1)); % init return var
    
    for w = 1:size(s.W,1) % loop over W values
        for i = 1:s.n     % loop over k values
            s.in(:,i,w) = loc_k(:,i) .* s.W(w,:)';
        end
    end
    
    % note s.corr shaped like model output, not inputs (s.in)
    s.corr = repmat([0;ones(floor(s.k/2),1)],1,s.n)';
    s.eps = loc;
 
end