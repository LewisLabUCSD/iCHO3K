function [points, samplingTime] = ADSB_v1(sample, verbose)
% Adaptive Directions Sampling on a Box (ADSB_v1)
%
% Uses ADSB to generate a uniform random sample from the loopless flux
% solution space with simplified boundary checking and robust loopless flux validation.
%
% INPUTS:
%   sample (structure): Structure containing sampling parameters
%   verbose (boolean): Controls the verbosity of the function
%
% OUTPUT:
%   points: Matrix with numRxns x numSamples loopless flux solutions
%   samplingTime: Time taken for the sampling process

    % Check input arguments
    if nargin < 2
        verbose = 1;
    end

    % Define helper functions
    mappingFxn = @(x_0, udir, step) x_0 + bsxfun(@times, udir, step);
    stepSizeFxn = @(x_0, L, R) x_0 + L .* R;

    % Initialize variables
    points = sample.points;
    nDim = size(points, 2);
    udir = zeros(size(points));
    currPoint = udir;
    nextPoint = currPoint;

    % Calculate number of steps
    nTimes = sample.stepsPerPoint;
    ptarget = 0.99;
    nSteps = getNumberSteps(nTimes, nDim, ptarget);
    nSteps = round((0:0.1:1) * nSteps);
    indexes = repmat(1:nDim, sample.numChains, 1, 1);

    % Pre-allocate memory
    posStep = zeros(1, sample.numChains);
    negStep = zeros(1, sample.numChains);
    Lcord = zeros(1, sample.numChains);
    steps = zeros(1, sample.numChains);

    % Start clock and sampling
    iterSample = 0;
    t0 = cputime;

    if verbose
        fprintf('---------------------------\n%%Prog \t Time \t Time left\n---------------------------\n');
    end

    while iterSample <= nSteps(end)
        % Update iteration counter
        iterSample = iterSample + 1;

        % Print progress
        if verbose && any(~(nSteps - iterSample))
            timeElapsed = (cputime - t0) / 60;
            timePerStep = timeElapsed / iterSample;
            fprintf('%d%%\t%8.2f min\t%8.2f min left\n', round(100 * iterSample / nSteps(end)), timeElapsed, (nSteps(end) - iterSample) * timePerStep);
        end

        % Sample new direction
        for ix = find(true(1, sample.numChains))  % Simplified check

            % Update current point and direction
            currPoint(:, ix) = points(:, indexes(ix, end), ix);
            udir(:, ix) = points(:, indexes(ix, end-2), ix) - points(:, indexes(ix, end-1), ix);
            if all(udir(:, ix) == 0), continue; end
            udir(:, ix) = udir(:, ix) / norm(udir(:, ix));

            % Determine step sizes to boundaries
            posStepTemp = (sample.ub - currPoint(:, ix)) ./ udir(:, ix);
            negStepTemp = (sample.lb - currPoint(:, ix)) ./ udir(:, ix);

            posStep(ix) = min([posStepTemp(udir(:, ix) > sample.uTol); negStepTemp(udir(:, ix) < -sample.uTol)]);
            negStep(ix) = max([negStepTemp(udir(:, ix) > sample.uTol); posStepTemp(udir(:, ix) < -sample.uTol)]);

            Lcord(ix) = posStep(ix) - negStep(ix);

            if Lcord(ix) < sample.bTol || posStep(ix) < 0 || negStep(ix) > 0
                continue;  % Skip invalid step
            end

            % Sample a random step and update the point
            steps(ix) = stepSizeFxn(negStep(ix), Lcord(ix), rand());
            nextPoint(:, ix) = mappingFxn(currPoint(:, ix), udir(:, ix), steps(ix));

            % Ensure points stay within bounds
            nextPoint(:, ix) = sample.keepWithinBounds(nextPoint(:, ix));

            % Check for loopless feasibility
            if sample.loopless
                if ~sample.isFeasible(nextPoint(:, ix))
                    continue;  % Skip non-loopless points
                end
            end

            % Update the point in the sample matrix
            points(:, indexes(ix, end), ix) = nextPoint(:, ix);
        end
    end

    % Return the final sampled points and sampling time
    samplingTime = (cputime - t0) / 60;
end
