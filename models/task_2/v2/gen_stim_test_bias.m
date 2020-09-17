function [stim] = gen_stim_test_bias(stim)
    % Generates diagonal Stimulus (e_t) at ones
    % size(stim) = (k: frames, 1, W: corresponding order)
    
    % generate W vec
    W = ones(stim.k,1,2^stim.k);
    for w = 1:2^stim.k
        c = dec2bin(w-1,stim.k);
        for i=1:stim.k
            W(i,1,w) = W(i,1,w) * sign(str2num(c(i))-0.5);
        end
    end
    
    [~, i] = sort(-abs(sum(W,1)));
    W = W(:,:,i);
    
    stim.eps = [-1,1];
    stim.W = W;
    stim.in = W;
    stim.corr = [ones(2^stim.k-2,1);0;0];
end