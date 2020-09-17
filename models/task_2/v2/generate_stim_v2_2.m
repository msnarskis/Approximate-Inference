%generate_stimulus_v2([0, 5], 10, 3, [1 1 1; 1 -1 1; -1 -1 1; -1 -1 -1])
function [s] = generate_stim_v2_2(s)
    % Generates (not noisy) Stimulus (e_t) for 
    % all values of W over n locations
    % size(stim) = (k: frames, n: location, W: corresponding order)
    
    % generate W vec
    s.W = ones(2^s.k,s.k);
    for w = 1:2^s.k
        c = dec2bin(w-1,s.k);
        for i=1:s.k
            s.W(w,i) = s.W(w,i) * sign(str2num(c(i))-0.5);
        end
    end
    
    [~, i] = sort(abs(sum(s.W,2)));
    s.W = s.W(i,:);    
    
    loc = linspace(s.eps_range(1), s.eps_range(2), s.n); % abs val of locations
    loc_k = repmat(loc, s.k, 1); % repeated
    
    s.in = zeros(s.k, s.n, size(s.W, 1)); % init return var
    
    for w = 1:size(s.W,1) % loop over W values
        for i = 1:s.n     % loop over k values
            s.in(:,i,w) = loc_k(:,i) .* s.W(w,:)';
        end
    end
    
    s.corr = repmat([ones(2^s.k-2,1);0;0],1,s.n)';
    s.eps = loc;
 
end