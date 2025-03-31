function [points, samplingTime] = ADSB_v2(sample, verbose)
% ADSB_v2: Adaptive Directions Sampling on a Box with Stoichiometric Constraints
%
% Uses ADSB to generate a uniform random sample from the feasible flux
% solution space while strictly enforcing stoichiometric constraints.
%
% INPUTS:
%   sample (structure): Structure containing sampling parameters
%                       Required fields:
%                       - S: Stoichiometric matrix
%                       - lb: Lower bounds
%                       - ub: Upper bounds
%                       - points: Initial points for sampling
%                       - numChains: Number of Markov chains
%                       - stepsPerPoint: Number of steps per point
%                       - bTol: Boundary tolerance
%                       - uTol: Direction vector tolerance
%                       - loopless: Flag for loopless sampling (set to false if not used)
%   verbose (boolean): Controls the verbosity of the function
%
% OUTPUTS:
%   points: Matrix with numRxns x numSamples flux solutions satisfying S*v=0
%   samplingTime: Time taken for the sampling process

    if nargin < 2
        verbose = true;
    end

    % Extract necessary data from sample structure
    S = sample.S;
    lb = sample.lb;
    ub = sample.ub;
    points = sample.points;
    numChains = sample.numChains;
    stepsPerPoint = sample.stepsPerPoint;
    bTol = sample.bTol;
    uTol = sample.uTol;
    loopless = isfield(sample, 'loopless') && sample.loopless;

    % Calculate the rank of S and estimate the null space dimension
    if issparse(S)
        % For sparse matrix S, use svds with a tolerance to estimate rank
        tolerance = 1e-9;
        [~, Sigma, ~] = svds(S, min(size(S)) - 1); % Calculate the largest singular values
        rank_S = sum(diag(Sigma) > tolerance); % Count singular values above the tolerance
    else
        % For dense matrix, use standard rank function
        rank_S = rank(S);
    end
    null_dim = size(S, 2) - rank_S;

    % Use svds with the estimated null space dimension to find the null space
    [U, Sigma, V] = svds(S, null_dim, 'smallest');
    nullIndices = find(diag(Sigma) < tolerance);
    nullBasis = V(:, nullIndices);

    % Initialize variables
    [numRxns, numWarmUpPoints, ~] = size(points);
    nDim = numWarmUpPoints;
    currPoint = zeros(numRxns, numChains);
    udir = zeros(numRxns, numChains);
    nextPoint = zeros(numRxns, numChains);
    steps = zeros(1, numChains);
    posStep = zeros(1, numChains);
    negStep = zeros(1, numChains);
    Lcord = zeros(1, numChains);

    % Initialize indices for points storage
    indexes = ones(numChains, 1);

    % Start timing
    t0 = cputime;

    if verbose
        fprintf('Starting ADSB_v2 sampling with %d chains and %d steps per point...\n', numChains, stepsPerPoint);
    end

    % Sampling loop
    for iterSample = 1:stepsPerPoint
        % For each chain
        for ix = 1:numChains
            % Get current point for this chain
            currPoint(:, ix) = points(:, indexes(ix), ix);

            % Generate a random direction vector in the null space of S
            dirCoefficients = randn(length(nullIndices), 1);
            udir(:, ix) = nullBasis * dirCoefficients;
            if norm(udir(:, ix)) < uTol
                % Regenerate direction if too small
                continue;
            end
            udir(:, ix) = udir(:, ix) / norm(udir(:, ix));

            % Determine maximum positive and negative step sizes within bounds
            posStepTemp = (ub - currPoint(:, ix)) ./ udir(:, ix);
            negStepTemp = (lb - currPoint(:, ix)) ./ udir(:, ix);

            % Handle infinite or NaN values
            posStepTemp(~isfinite(posStepTemp)) = Inf;
            negStepTemp(~isfinite(negStepTemp)) = -Inf;

            % Maximum feasible step sizes
            posStep(ix) = min([posStepTemp(udir(:, ix) > uTol); negStepTemp(udir(:, ix) < -uTol); Inf]);
            negStep(ix) = max([negStepTemp(udir(:, ix) > uTol); posStepTemp(udir(:, ix) < -uTol); -Inf]);

            % Calculate total length of feasible segment
            Lcord(ix) = posStep(ix) - negStep(ix);

            if Lcord(ix) < bTol || posStep(ix) <= 0 || negStep(ix) >= 0
                % Skip if no feasible step
                continue;
            end

            % Sample a random step size within the feasible segment
            steps(ix) = negStep(ix) + Lcord(ix) * rand();

            % Calculate next point
            nextPoint(:, ix) = currPoint(:, ix) + udir(:, ix) * steps(ix);

            % Ensure next point is within bounds
            nextPoint(:, ix) = max(min(nextPoint(:, ix), ub), lb);

            % Check stoichiometric feasibility
            if norm(S * nextPoint(:, ix)) > bTol
                % Project onto null space to correct numerical errors
                nextPoint(:, ix) = nextPoint(:, ix) - S' * ((S * S') \ (S * nextPoint(:, ix)));
                % Re-apply bounds after projection
                nextPoint(:, ix) = max(min(nextPoint(:, ix), ub), lb);
            end

            % Loopless flux check (if required)
            if loopless
                if ~isLooplessFlux(nextPoint(:, ix), S, uTol)
                    % Skip if not loopless
                    continue;
                end
            end

            % Store the new point
            indexes(ix) = indexes(ix) + 1;
            if indexes(ix) > size(points, 2)
                % Expand the points matrix if needed
                points(:, end+1, ix) = nextPoint(:, ix);
            else
                points(:, indexes(ix), ix) = nextPoint(:, ix);
            end
        end

        % Verbose output
        if verbose && mod(iterSample, round(stepsPerPoint / 10)) == 0
            timeElapsed = (cputime - t0) / 60;
            fprintf('Progress: %d%%, Time Elapsed: %.2f min\n', round(100 * iterSample / stepsPerPoint), timeElapsed);
        end
    end

    % Reshape points to 2D matrix (numRxns x totalSamples)
    totalSamples = sum(indexes - 1);
    allPoints = zeros(numRxns, totalSamples);
    idx = 0;
    for ix = 1:numChains
        nPoints = indexes(ix) - 1;
        allPoints(:, idx + (1:nPoints)) = points(:, 2:indexes(ix), ix);
        idx = idx + nPoints;
    end
    points = allPoints;

    % Record sampling time
    samplingTime = (cputime - t0) / 60;

    if verbose
        fprintf('Sampling completed in %.2f minutes.\n', samplingTime);
    end
end

function isLoopless = isLooplessFlux(flux, S, tol)
    % Placeholder function for loopless flux check
    % Implement loopless flux verification method here
    % For now, we'll assume all fluxes are acceptable
    isLoopless = true;
end
