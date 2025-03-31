function nSteps = getNumberSteps(nTimes,nDim,ptarget)
% Estimate number of moves for ADSB
%
% Determine number of steps in order to move each particle nTimes with
% at least ptarget probability (e.g. 99%) 
%
% USAGE:
%              nSteps = getNumberSteps(nTimes,nDim,ptarget)
%
% INPUTS:
%              nTimes - Number of times we want to move a particle (integer)
%              nDim   - Population size (integer)
%              ptarget - Target probability
%
% OUTPUT:
%              nSteps    Number of Matrix with numRxns x numSamples (loopless) flux solutions
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

p_mov      = 1/nDim;                                                       % Probability of choosing a member from the population
nSteps     = 5*nTimes;                                                     % Initial nSteps guess
logp_mov   = log(p_mov);
logp_nmov  = log(1-p_mov);
log_kprob  = zeros(nTimes,1);
log_cumsum = cumsum(log(1:nSteps));

% Main recursive loop
while true   
    log_kprob(1) = nSteps*logp_nmov;                                       % Deal with special case k = 0
    for k = 1:nTimes-1                                                     % Calculate Binomial probs up to (nTimes-1) successes
        log_kprob(k+1) = log_cumsum(nSteps) - log_cumsum(nSteps-k) - log_cumsum(k) + k*logp_mov + (nSteps-k)*logp_nmov;
    end
    max_logkp  = max(log_kprob);                                           % Avoid underflow
    prob_k     = exp(max_logkp)*sum(exp(log_kprob - max_logkp));
    prob_kplus = 1 - prob_k;
    if (prob_kplus<ptarget)
        nSteps = nSteps+1;
        log_cumsum = [log_cumsum,log_cumsum(end)+log(nSteps)];
    else
        break;
    end
end
