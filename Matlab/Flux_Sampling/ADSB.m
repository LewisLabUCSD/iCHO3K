function [points,samplingTime] = ADSB(sample,verbose)
% Adaptive Directions Sampling on a Box
%
% Uses ADSB to generate a uniform random sample from the loopless flux
% solution space
%
% USAGE:
%              points = ADSB(sample)
%
% INPUTS:
%              sample (structure):    (the following fields are required - others can be supplied)
%                                    * pointsPerChain - Number of points for the chains (double)
%                                    * stepsPerPoint - Thining or number of steps per effective point (double)
%                                    * numChains - Number of Markov chains (double)
%                                    * numRxns - Number of rxns
%                                    * numSamples - Total number of samples
%                                    * loopless - Loop removal option, true (default) or false
%                                    * stepFxn - Step function option, 'uniform' (default) or 'linear'
%                                    * x0 - `n x nDim` Matrix of initial points of the current set
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * keepWithinBounds - Fxn handle for bringToBoundary.m fxn
%                                    * isFeasible - Fxn handle for looplessCheck.m fxn
%
% OUTPUT:
%              points:   Matrix with numRxns x numSamples (loopless) flux solutions
%
% OPTIONAL OUTPUT:
%              samplingTime:   Runtime of ADSB
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Check input arguments
if nargin<2; verbose = 1; end

% Define mapping fxn for new points
mappingFxn = @(x_0,udir,step) x_0 + bsxfun(@times,udir,step);

% Build step-calculating function
stepSizeFxn = @(x_0,L,R) x_0 + L.*R;

% Initialize sampling variables
points    = sample.points;
nDim      = size(points,2);
udir      = zeros(size(points));
currPoint = udir;
nextPoint = currPoint;

% Figure out number of steps and pre-allocate
nTimes    = sample.stepsPerPoint;                                                     % Number of times we expect to move a particle (e.g. 100)
ptarget   = .99;                                                                      % Probability (target) of moving a member at least nTimes (e.g. 99%)
nSteps    = getNumberSteps(nTimes,nDim,ptarget);
nSteps    = round((0:.1:1)*nSteps);
indexes   = repmat(1:nDim,sample.numChains,1,1);

% Pre-allocate memory
Lcord   = zeros(1,sample.numChains);
steps   = zeros(1,sample.numChains);
posStep = zeros(1,sample.numChains);
negStep = zeros(1,sample.numChains);

% Start clock and sampler
iterSample = 0;
t0         = cputime;                                                                   % Start clock
if (verbose==1)
    fprintf('---------------------------\n%%Prog \t Time \t Time left\n---------------------------\n');
end
while (iterSample <= nSteps(end))
    
    % Update general counter
    iterSample = iterSample+1;
    
    % Print progress information
    if (iterSample>1) && any(~(nSteps-iterSample)) && (verbose==1)
        timeElapsed = (cputime-t0)/60;
        timePerStep = timeElapsed/iterSample;
        fprintf('%d\t%8.2f\t%8.2f\n',round(1e2*iterSample/nSteps(end)),timeElapsed,(nSteps(end)-iterSample)*timePerStep);
    end
    
    % Pre-allocate variables
    flags = true(1,sample.numChains);
        
    % Sample feasible direction
    while any(flags)       
        idx_next = indexes;
        
        % Figure out the true max & min step sizes
        for ix = find(flags)
            
            % Update points based on the previous iteration
            if (iterSample>1); points(:,idx_prev(ix,end),ix) = nextPoint(:,ix); end
            
            % Sample three indexes without replacement
            idx1 = randi(nDim);
            idx_next(ix,[idx1,nDim])   = idx_next(ix,[nDim,idx1]);
            idx2 = randi(nDim-1);
            idx_next(ix,[idx2,nDim-1]) = idx_next(ix,[nDim-1,idx2]);
            idx3 = randi(nDim-2);
            idx_next(ix,[idx3,nDim-2]) = idx_next(ix,[nDim-2,idx3]);          
            
            % Define current point and direction of movement
            currPoint(:,ix) = points(:,idx_next(ix,end),ix);
            udir(:,ix)      = points(:,idx_next(ix,end-2),ix)-points(:,idx_next(ix,end-1),ix);
            if all(udir(:,ix)==0); continue; end
            udir(:,ix)      = udir(:,ix)/norm(udir(:,ix));                                  % unitary vector
            
            % Figure out positive and negative directions, and maximum distance to boundaries
            posDir      = (udir(:,ix)>sample.uTol);
            negDir      = (udir(:,ix)<-sample.uTol);
            posStepTemp = (sample.ub-currPoint(:,ix))./udir(:,ix);
            negStepTemp = (sample.lb-currPoint(:,ix))./udir(:,ix);
            posStep(ix) = min([posStepTemp(posDir);negStepTemp(negDir)]);
            negStep(ix) = max([negStepTemp(posDir);posStepTemp(negDir)]);
        end
        
        % Draw new direction if too close to the boundaries
        Lcord(flags) = posStep(flags)-negStep(flags);
        flags(flags) = (Lcord(flags) < sample.bTol) | (posStep(flags) < 0) | (negStep(flags) > 0);
    end
    
    % Sample from the shrinking interval as in slice sampling
    flags = true(1,sample.numChains);
    while any(flags)
        
        % Perform a random step
        steps(flags)       = stepSizeFxn(negStep(flags),Lcord(flags),rand(1,sum(flags)));
        nextPoint(:,flags) = mappingFxn(currPoint(:,flags),udir(:,flags),steps(flags));
        nextPoint(:,flags) = sample.keepWithinBounds(nextPoint(:,flags));
        if sample.loopless
            flags(flags)   = ~sample.isFeasible(nextPoint(:,flags));
            
            % Shrink hypercube if necessary
            posSign        = (steps(flags)>0);
            posStep(flags) = posStep(flags).*(~posSign) + steps(flags).*(posSign);
            negStep(flags) = negStep(flags).*(posSign)  + steps(flags).*(~posSign);
            
            % Keep shrinking only if the bounding box is not too small
            Lcord(flags) = posStep(flags)-negStep(flags);
            flags(flags) = (Lcord(flags)>sample.bTol);
        else
            break;
        end
    end
    
    % Reset position of the next point to the current point for zero-length steps
    fixed_point = (Lcord<sample.bTol);
    nextPoint(:,fixed_point) = currPoint(:,fixed_point);
    
    % Save indexes to change in the next iteration
    idx_prev = idx_next;
end
samplingTime = (cputime-t0)/60;