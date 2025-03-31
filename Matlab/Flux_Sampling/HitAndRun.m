function [points,samplingTime,rejectionRate] = HitAndRun(sample,numSamples,numDiscarted,stepsPerPoint,verbose,workerIdx,initPoint)
% Hit and Run algorithm
%
% Implements HR to generate a uniform random flux sample. The loopless
% condtion can be enforced by rejecting loopy samples on-the-fly
%
% USAGE:
%              points = HitAndRun(sample,numSamples,numDiscarted,stepsPerPoint)
%
% INPUTS:
%              sample (structure):    (the following fields are required - others can be supplied)
%                                    * COBRA model structure fields (S, ub, lb)
%                                    * centroid - Initial centroid estimation
%                                    * isFeasible - Fxn handle for looplessCheck.m fxn
%                                    * loopless - Loop removal option, true (default) or false
%              numSamples:           Number of points to sample
%              numDiscarded:         Burn-in (double)
%              stepsPerPoint:        Thinning (double)
%
% OPTIONAL INPUTS:
%              verbose:    Display progress of sampling run, true (default) or false
%              workerIdx:  ID of the working processor, 1 (default)
%              initPoint:  Initial point, centroid (default)
%
% OUTPUT:
%              points:   Matrix with n x numSamples (loopless) flux solutions
%
% OPTIONAL OUTPUT:
%              samplingTime:   Runtime of HRB
%              rejectionRate:  Rejection rate estimated by the number of
%                              times the 'hyperbox' is shrunk
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Check inputs
if nargin<4; disp('Not enough input parameters'); return; end

% Set up Hit-And-Run parameters
t0 = cputime;
if nargin<5
    verbose   = true;
    workerIdx = 1;
    initPoint = sample.centroid;
elseif nargin<6
    workerIdx = 1;
    initPoint = sample.centroid;
elseif nargin<7
    initPoint = sample.centroid;
end

% Allocating memory for the sampling points
nRxns  = size(sample.S,2);
points = zeros(nRxns,numSamples);
% N      = fast_snp(sample.S,sample.lb,sample.ub,sample.vTol);
N      = null(sample.S,'r');
N      = bsxfun(@rdivide,N,sqrt(sum(N.^2,1)));
fVars  = size(N,2);

% Define sampler parameters
uTol      = sample.uTol;               % Direction tolerance
maxMinTol = sample.bTol;               % Min distance to closest bound
% dTol      = sample.bTol;               % Max discrepancy with equality constraints

% Set initial point
prevPoint = initPoint;

% Initialize counters
counter       = 0;
sampleCount   = 0;
totalCount    = 0;
rejectedCount = 0;
rejectionRate = 0;

% Main loop
if (verbose==1) && (workerIdx==1)
    countReport  = 1;
    sampleReport = 0;
    fprintf('--------------------------------------------\n%%Prog \t Time \t Time left \t Rejection rate\n--------------------------------------------\n');
end

while sampleCount<numSamples

    % Update count
    totalCount    = totalCount+1;
    rejectionRate = rejectedCount/totalCount;

    % Print step information
    if (verbose==1) && (workerIdx==1) && (sampleReport==sampleCount)
        timeElapsed  = (cputime-t0)/60;
        timePerStep  = timeElapsed/sampleCount;
        if  (counter>0)
            fprintf('%d\t%8.2f\t%8.2f\t%8.2f\n',round(100*sampleReport/numSamples),timeElapsed,(numSamples-sampleCount)*timePerStep,rejectionRate);
        elseif (counter==0)
            fprintf('%d\t%8.2f\t%8.2f\t%8.2f\n',round(100*sampleReport/numSamples),timeElapsed,(numSamples-sampleCount)*timePerStep,1);
        end
        countReport  = countReport+1;
        sampleReport = min([numSamples,countReport*round(numSamples/20)]);
    end

    % Return if the dynamics is frozen
    if rejectionRate>.9999 && totalCount>1e5
        disp('Hit-And-Run dynamics is frozen!');
        points = []; break;
    end

    % Sample random direction (Coordinate Kernel Hit-And-Run)
    u = N(:,randi(fVars));

    % Figure out the distances to upper and lower bounds
    distUb = (sample.ub-prevPoint);
    distLb = (prevPoint-sample.lb);

    % Figure out positive and negative directions
    posDirn = (u>uTol);
    negDirn = (u<-uTol);

    % Figure out all the possible maximum and minimum step sizes
    maxStepTemp = distUb./u;
    minStepTemp = -distLb./u;
    maxStepVec  = [maxStepTemp(posDirn);minStepTemp(negDirn)];
    minStepVec  = [minStepTemp(posDirn);maxStepTemp(negDirn)];

    % Figure out the true max & min step sizes
    maxStep = min(maxStepVec);
    minStep = max(minStepVec);

    % Find new direction if we're getting too close to a constraint
    if (abs(minStep)<maxMinTol && abs(maxStep)<maxMinTol) || (minStep>maxStep) || (maxStep<0) || (minStep>0)
        rejectedCount = rejectedCount+1;
        continue;
    end

    % Pick a random step distance
    stepDist = rand(1)*(maxStep-minStep)+minStep;

    % Make move and check whether it is feasible
    nextPoint = prevPoint + stepDist*u;

    % Bring points to the boundary if they fall out
    if any((sample.ub(1:nRxns)-nextPoint)<0)
        overInd = (nextPoint>sample.ub(1:nRxns));
        nextPoint(overInd) = sample.ub(overInd);
    end
    if any((nextPoint-sample.lb(1:nRxns))<0)
        underInd = (nextPoint<sample.lb(1:nRxns));
        nextPoint(underInd) = sample.lb(underInd);
    end

    % Check loopless condition topologically
    if sample.loopless
        looplessModel = 0;
        if sample.isFeasible(nextPoint)
            looplessModel = 1;
        end
    else
        looplessModel = 1;
    end

    % Save current point if it is loopless
    if looplessModel
        counter   = counter+1;
        prevPoint = nextPoint;

        % Save sampled point
        if counter>numDiscarted && ~mod(counter,stepsPerPoint)
            sampleCount = sampleCount+1;
            points(:,sampleCount) = nextPoint;
        end
    else
        % Update rejected count
        rejectedCount = rejectedCount+1;

    end
end
samplingTime = (cputime-t0)/60;