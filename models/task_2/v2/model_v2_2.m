function [resp] = model_v2_2(stim, cond, sig_t, sig_n, sig_v, sig_sa, sig_sv, pr_R, pr_C, noise_a, noise_v, nsamp, ntrial)
% run task_2 model, not made to run with small RAM
% returns matrix w/ probability of answering R=1
%   resp size = (n: stim locations, w: experimental W's seperated for analysis)
%   matched:  cond = 1 
%   center :  cond = 0

    % size params for ease
    k = size(stim,1);
    n = size(stim,2);
    w = size(stim,3);    
    
    % generate W vec for inference
    W = ones(k,2^k);
    for w = 1:2^k
        c = dec2bin(k,w-1);
        for i=1:k
            W(i,w) = W(i,w) * sign(str2num(c(i))-0.5);
        end
    end
    
    % generate R vec for inference
    R = [1 zeros(1,2^k-2) 1];    
    
    % init vars
    resp = zeros(n,size(stim,3));
    zero = zeros(k,n);
    
    % repeat for trials
    X = repmat(stim,1,1,1,ntrials); % size=(k,n,w,ntrials)
    
    % add noise
    X_t =    stim(:,:,w) + randn(k,n,w,ntrials)*noise_a;
    X_n = -1*stim(:,:,w) + randn(k,n,w,ntrials)*noise_a;
    X_R =    abs(stim(:,:,w))*cond + randn(k,n,w,ntrials)*noise_v;
    X_L = -1*abs(stim(:,:,w))*cond + randn(k,n,w,ntrials)*noise_v;
    
    % compute closed-form variables
    a_tn = sig_n/(sig_t + sig_n);
    X_tn = X_t.*a_tn - X_n.*(1-a_tn);
    sig_tn = sig_t * a_tn;
            
    a_stn = sig_sa/(sig_tn + sig_sa);
    X_stn = X_tn*a_stn;
    sig_stn = sig_tn * a_stn;
            
    a_RL = 1/2;
    X_RL = X_R.*a_RL - X_L.*(1-a_RL);
    sig_RL = sig_v * a_RL;
            
    a_sRL = sig_sv/(sig_RL + sig_sv);
    X_sRL = X_RL*a_sRL;
    sig_sRL = sig_RL * a_sRL;
    
    p_R0 = zeros(size(W,1),k,n,w,ntrials);
    p_R1 = zeros(size(W,1),k,n,w,ntrials);
    
    % marginalize over Ws
    for i = 1:size(W,1)
        Wi = repmat(W(:,i),1,n,w,ntrials);
        
        % Wi dependant variables
        a_tnRL = sig_sRL/(sig_stn + sig_sRL);
        X_tnRL = X_stn.*a_tnRL + Wi.*X_sRL.*(1-a_tnRL);
        sig_tnRL = sig_stn * a_tnRL;
                
        % normalizing stuff
        gamn = (1-pr_C)*0.25;
        gam = (pr_C)*normpdf(0,0,sqrt(sig_sa + sig_sv))*0.5;
        
        Z =   lognormpdf(X_tn, zero, sqrt(sig_tn + sig_sa))...
            + lognormpdf(X_RL, zero, sqrt(sig_RL + sig_sv))...
            + lognormpdf(X_t, -X_n,  sqrt(sig_t  + sig_n ))...
            + lognormpdf(X_R, -X_L,  sqrt(sig_v  + sig_v ));
        
        % no combine term
        tnC = log(normcdf(zero, -Wi.*X_stn, sqrt(sig_stn)))...
              + log(normcdf(zero, -X_sRL, sqrt(sig_sRL)))...
              + Z + log(1-pr_C) - log(gam + gamn);
                
        % combine term
        tC =  log(normcdf(zero, -Wi.*X_tnRL, sqrt(sig_tnRL)))...
              + lognormpdf(X_stn, Wi.*X_sRL, sqrt(sig_stn + sig_sRL))...
              + Z + log(pr_C) - log(gam + gamn);
                
        tmp(:,:,:,:,1)=tnC;
        tmp(:,:,:,:,2)=tC;
        
        if R(i) % W -> R deterministic
            p_R1(i,:,:,:,:) = sum(logsumexp(tmp,5),1)-log(2);
        else
            p_R0(i,:,:,:,:) = sum(logsumexp(tmp,5),1)-log(2^k-2);
        end
        
    end
    
    
    

end
