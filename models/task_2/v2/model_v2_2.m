function [resp] = model_v2_2(stim, cond, sig_t, sig_n, sig_v, sig_sa, sig_sv, pr_R, pr_C, noise_a, noise_v, nsamp, ntrials)
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
    for j = 1:2^k
        c = dec2bin(j-1,k);
        for i=1:k
            W(i,j) = W(i,j) * sign(str2num(c(i))-0.5);
        end
    end
    
    % generate R vec for inference
    R = [1 zeros(1,2^k-2) 1];    
    
    % init vars
    resp = zeros(n,size(stim,3));
    zero = zeros(k,n,w,ntrials);
    
    % repeat for trials
    X = repmat(stim,1,1,1,ntrials); % size=(k,n,w,ntrials)
    
    % add noise
    X_t =    X + randn(k,n,w,ntrials)*noise_a;
    X_n = -1*X + randn(k,n,w,ntrials)*noise_a;
    X_R =    abs(X)*cond + randn(k,n,w,ntrials)*noise_v;
    X_L = -1*abs(X)*cond + randn(k,n,w,ntrials)*noise_v;
    
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
    
    % choice vars
    l_R0 = zeros(size(W,2),n,w,ntrials);
    l_R1 = zeros(size(W,2),n,w,ntrials);
    
    % marginalize over W
    for i = 1:size(W,2)
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
        
        % size = (size(W,2),n,w,ntrials,2)
        tmp(:,:,:,:,1)=tnC;
        tmp(:,:,:,:,2)=tC;
        
        % product across frames
        if R(i) % W -> R deterministic
            l_R1(i,:,:,:) = sum(logsumexp(tmp,5),1)-log(2);
        else
            l_R0(i,:,:,:) = sum(logsumexp(tmp,5),1)-log(2^k-2);
        end
        
    end % end marginalize over W
    
    % get rid of unfilled rows
    % note: p_RX size = (size(W,1),n,w,ntrials)
    l_R0 = l_R0(2:end-1,:,:,:);
    l_R1 = l_R1([1,end],:,:,:);
    
    % combine w prior (1,n,w,ntrials)
    p_R0 = logsumexp(l_R0,1) + log(1-pr_R);
    p_R1 = logsumexp(l_R1,1) + log(pr_R);
    
    % calculate denomenator
    tmp = zeros(2,n,w,ntrials);
    tmp(1,:,:,:,:) = p_R1;
    tmp(2,:,:,:,:) = p_R0;
    
    den = logsumexp(tmp,1);
    
    % choice variable
    p_R = exp(p_R1 - den);
    
    % sample and return
    p_Rs = sampling(p_R, nsamp);
    
    resp = squeeze(mean(p_Rs,4));

end

function [resp] = sampling(inp, nsamp)
        
    if any(isnan(inp),'all') || any(isinf(inp),'all')
        ["num error:", inp]
    end
    
    inp(isnan(inp)) = .5;
    
    % nsamples < Inf
    if nsamp~=Inf
        resp = binocdf(floor(nsamp/2), nsamp, inp, 'upper')...
        + .5*binopdf(nsamp/2, nsamp, inp);
    end

    % nsamples = Inf
    if nsamp==Inf
        resp = (inp > 0.5) + 0.5*(inp == 0.5);
    end
end
