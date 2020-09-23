classdef Model_Task_2
    properties % Parameters
        % noise
        var_t; 
        var_n;
        var_v;
        var_sa;
        var_sv;
        
        % priors
        pr_R;
        pr_C;
        
        nsamp; % approx degree
        
        % simulation parameters
        noisy_in;
        ntrials;
    end
    
    methods
        % constructor
        function obj = Model_Task_2(par)
             % noise
            obj.var_t  = par.var_t; 
            obj.var_n  = par.var_n;
            obj.var_v  = par.var_v;
            obj.var_sa = par.var_sa;
            obj.var_sv = par.var_sv;

            % priors
            obj.pr_R  = par.pr_R;
            obj.pr_C  = par.pr_C;

            obj.nsamp = par.nsamp; % approx degree

            % simulation parameters
            obj.noisy_in = par.noisy_in;
            obj.ntrials  = par.ntrials*par.noisy_in + (1-par.noisy_in);
        end
        
        function [resp] = simulate(obj,stim,cond)
        % run task_2 model, not made to run with small RAM
        % returns matrix w/ probability of answering R=1
        %   resp size = (n: stim locations, w: experimental W's seperated for analysis)
        %   matched:  cond = 1 
        %   center :  cond = 0
        %
        % stim, cond, var_t, var_n, var_v, var_sa, var_sv, pr_R, pr_C, ...
        % noise_a, noise_v, nsamp, ntrials

            % size params for ease
            k = size(stim.in,1);
            n = size(stim.in,2);
            w = size(stim.in,3);    

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
            resp = zeros(n,size(stim.in,3));
            zero = zeros(k,n,w,obj.ntrials);

            % repeat for trials
            X = repmat(stim.in,1,1,1,obj.ntrials); % size=(k,n,w,ntrials)

            % add noise
            X_t =    X + randn(k,n,w,obj.ntrials)*sqrt(obj.var_t)*obj.noisy_in;
            X_n = -1*X + randn(k,n,w,obj.ntrials)*sqrt(obj.var_n)*obj.noisy_in;
            X_R =    abs(X)*cond + randn(k,n,w,obj.ntrials)*sqrt(obj.var_v)*obj.noisy_in;
            X_L = -1*abs(X)*cond + randn(k,n,w,obj.ntrials)*sqrt(obj.var_v)*obj.noisy_in;

            % compute closed-form variables
            a_tn = obj.var_n/(obj.var_t + obj.var_n);
            X_tn = X_t.*a_tn - X_n.*(1-a_tn);
            var_tn = obj.var_t * a_tn;

            a_stn = obj.var_sa/(var_tn + obj.var_sa);
            X_stn = X_tn*a_stn;
            var_stn = var_tn * a_stn;

            a_RL = 1/2;
            X_RL = X_R.*a_RL - X_L.*(1-a_RL);
            var_RL = obj.var_v * a_RL;

            a_sRL = obj.var_sv/(var_RL + obj.var_sv);
            X_sRL = X_RL*a_sRL;
            var_sRL = var_RL * a_sRL;

            % choice vars
            l_R0 = zeros(size(W,2),n,w,obj.ntrials);
            l_R1 = zeros(size(W,2),n,w,obj.ntrials);

            % marginalize over W
            for i = 1:size(W,2)
                Wi = repmat(W(:,i),1,n,w,obj.ntrials);

                % Wi dependant variables
                a_tnRL = var_sRL/(var_stn + var_sRL);
                X_tnRL = X_stn.*a_tnRL + Wi.*X_sRL.*(1-a_tnRL);
                sig_tnRL = var_stn * a_tnRL;

                % normalizing stuff
                % gamn = (1-pr_C)*0.25;
                % gam = (pr_C)*normpdf(0,0,sqrt(sig_stn + sig_sRL))*0.5;

                Z =   lognormpdf(X_tn, zero, sqrt(var_tn + obj.var_sa))...
                    + lognormpdf(X_RL, zero, sqrt(var_RL + obj.var_sv))...
                    + lognormpdf(X_t, -X_n,  sqrt(obj.var_t  + obj.var_n ))...
                    + lognormpdf(X_R, -X_L,  sqrt(obj.var_v  + obj.var_v ))...
                    ...
                    + lognormpdf(X, X_t, sqrt(obj.var_t))...
                    + lognormpdf(X, X_n, sqrt(obj.var_n))...
                    + lognormpdf(X*cond, X_R, sqrt(obj.var_v))...
                    + lognormpdf(X*cond, X_L, sqrt(obj.var_v));

                % no combine term
                tnC = log(normcdf(zero, -Wi.*X_stn, sqrt(var_stn)))...
                      + log(normcdf(zero, -X_sRL, sqrt(var_sRL)))...
                      + Z + log(1-obj.pr_C);% - log(gam + gamn);

                % combine term
                tC =  log(normcdf(zero, -Wi.*X_tnRL, sqrt(sig_tnRL)))...
                      + lognormpdf(X_stn, Wi.*X_sRL, sqrt(var_stn + var_sRL))...
                      + Z + log(obj.pr_C);% - log(gam + gamn);

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
            p_R0 = logsumexp(l_R0,1) + log(1-obj.pr_R);
            p_R1 = logsumexp(l_R1,1) + log(obj.pr_R);

            % calculate denomenator
            tmp = zeros(2,n,w,obj.ntrials);
            tmp(1,:,:,:,:) = p_R1;
            tmp(2,:,:,:,:) = p_R0;

            den = logsumexp(tmp,1);

            % choice variable
            p_R = exp(p_R1 - den);

            % sample and return
            p_Rs = obj.sampling(p_R, obj.nsamp);

            resp = squeeze(mean(p_Rs,4));
        end
        
        function [resp] = sampling(obj, inp, nsamp)
        
            if any(isnan(inp),'all') || any(isinf(inp),'all')
                ["num error:", sum(isnan(inp),'all')]
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
    end
end