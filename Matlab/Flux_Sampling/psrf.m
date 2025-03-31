function [R,Rint,Neff,tau,thin] = psrf(varargin)
% Potential Scale Reduction Factor based on Brooks and Gelman method.
% (adapted from Nils Winter's Matlab Toolbox for Bayesian Estimation. For
% further references and license information refer to 
% https://github.com/NilsWinter/matlab-bayesian-estimation)
%
% Performs Markov chain diagnostics and returns the PSRF for a collection
% of MCMC-simulations. If PSRF is not close to 1 (below 1.1 for example)
% one may conclude that the tested samples were not from the same
% distribution, and thus, the chain might not have been converged yet.
%
% USAGE:
%              [R,Rint,Neff] = psrf(X)
%            
% INPUTS:
%              X:  N x D x M matrix containing M MCMC simulations of
%                  length N, each with dimension D
%
% OUTPUT:
%              R:     PSRF (R=sqrt(V/W)) (1 x D vector)
%              Rint:  PSRF for the (10,90)% confidence interval region (1 x D vector)
%              Neff:  Estimated effective number of samples (1 x D vector)
%
% OPTIONAL OUTPUTS:
%              tau:   Estimated autocorrelation time (1 x D vector)
%              thin:  Geyer's initial positive sequence lag. Useful for
%                     estimating an `optimal` thinning factor (1 x D vector)
%
% -------------------- Adapted on 2023 by Pedro A. Saa --------------------

X   = cat(3,varargin{:});
mid = floor(size(X,1)/2);
X   = cat(3,X(1:mid,:,:),X((end-mid+1):end,:,:));
[N,D,M] = size(X);
if (N<=2); fprintf('Too few samples'); return; end

% Calculate means W of the variances
W = zeros(1,D);
for mi = 1:M
    x = bsxfun(@minus,X(:,:,mi),mean(X(:,:,mi)));
    W = W + sum(x.*x);
end
W = W / ((N-1) * M);

% Calculate variances B (in fact B/n) of the means.
Bpn = zeros(1,D);
m   = mean(reshape(mean(X),D,M)');
for mi =1:M
    x   = mean(X(:,:,mi)) - m;
    Bpn = Bpn + x.*x;
end
Bpn = Bpn / (M-1);

% Calculate reduction factors
B  = Bpn*N;
Vh = (N-1)/N*W + Bpn;
R  = sqrt(Vh./W);

% Calculate reduction factor based on 95% interval
if nargout>2
    W = zeros(1,D);
    V = zeros(1,D);
    for d = 1:D
        x = X(:,d,:);
        W(1,d) = mean(diff(prctile(x,[10,90])));
        V(1,d) = diff(prctile(x(:),[10,90]));
    end
    Rint = V./W;
end

% Compute autocorrelation
if nargout>3       
    
    %  Variogram (try to paralellize for speed up)    
    try
        try 
            parpool('local');
        catch
            delete(gcp('nocreate'));
            parpool('local');
        end
        parfor t = 1:N-1
            Vt(t,:) = sum(sum((X(1:end-t,:,:)-X(1+t:end,:,:)).^2,1),3)/M/(N-t);
        end
        delete(gcp('nocreate'));
    catch
        for t = 1:N-1
            Vt(t,:) = sum(sum((X(1:end-t,:,:)-X(1+t:end,:,:)).^2,1),3)/M/(N-t);
        end
    end
    
    % Autocorrelation
    rho  = 1-bsxfun(@rdivide,Vt./2,Vh);
    rho  = [ones(1,D);rho];                                                    % add zero lag autocorrelation
    mid  = floor(N/2);
    Neff = zeros(1,D);
    for di =1:D
        cp = sum(reshape(rho(1:2*mid,di),2,mid),1);
        ci = find(cp<0,1);
        if isempty(ci)                                                          % Inital positive could not be found for this variable, use maxlag value instead
            ci = mid;
        else
            ci = ci-1;                                                          % Last positive value found
        end
        cp       = [cp(1:ci) 0];                                                % Initial positive sequence
        tau(di)  = -1+2*sum(cp);                                                % Initial positive sequence estimator
        Neff(di) = M*N/tau(di);                                                 % Initial positive sequence estimator for neff
        thin(di) = ci*2;
    end
end
