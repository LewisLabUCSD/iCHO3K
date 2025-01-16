function sample = fluxSamplerADSB_v2(model, options)
% ADSB Flux Sampler
%
% Performs random sampling of the flux space of metabolic models using ADSB
%
% INPUTS:
%              model (structure):    (the following fields are required)
%                                    * S  - Stoichiometric matrix
%                                    * lb - Lower bounds
%                                    * ub - Upper bounds
%                                    * rxns - Reaction identifiers (cell array)
%              options (structure):  (the following fields are required)
%                                    * numSamples - Number of samples to be generated
%                                    * stepsPerPoint - Thinning or number of steps per effective point
%                                    * populationScale - Size of the set relative to the feasible space dimension
%
% OUTPUT:
%              sample (structure):   Sampling structure containing the random samples

    % Check if options are provided correctly
    if nargin < 2
        error('Not enough input arguments.');
    else
        if ~isfield(options, 'numSamples'); error('Field numSamples has not been defined.'); end
        if ~isfield(options, 'stepsPerPoint'); options.stepsPerPoint = 1e2; end
        if ~isfield(options, 'populationScale'); options.populationScale = 3; end
    end

    % Initialize sample structure with options
    sample.numSamples = options.numSamples;
    sample.stepsPerPoint = options.stepsPerPoint;
    sample.populationScale = options.populationScale;

    % Define solver parameters (assuming COBRA toolbox is being used)
    changeCobraSolverParams('LP', 'optTol', 1e-9);
    changeCobraSolverParams('LP', 'feasTol', 1e-9);

    % Load model information into sample structure
    sample.rxns = model.rxns;
    sample.S = model.S;
    sample.lb = model.lb;
    sample.ub = model.ub;

    % Set tolerance for flux direction
    sample.uTol = 1e-5;  % Define a small numerical tolerance value for flux direction
    sample.bTol = 1e-5;  % Define a small boundary tolerance value

    % Generate initial flux seeds using optimizeCbModel
    fprintf('Generating initial flux seeds using optimizeCbModel...\n');

    % Number of initial seeds to generate
    numSeeds = options.numSamples;

    % Initialize matrix to store initial flux seeds
    initialSeeds = zeros(length(model.rxns), numSeeds);

    % Run a single optimization to use as a basis for initial seeds
    sol = optimizeCbModel(model, 'max', 'one');
    
    if sol.stat == 1
        % Use the solution to generate multiple initial seeds with slight random perturbations
        baseSeed = sol.x;
        perturbationRange = 0.1; % Perturb each flux by up to ±10%
        for i = 1:numSeeds
            perturbation = (rand(size(baseSeed)) - 0.5) * 2 * perturbationRange .* baseSeed;
            initialSeeds(:, i) = baseSeed + perturbation;
        end
        sample.warmUpPoints = initialSeeds;
    else
        fprintf('Warning: Initial optimization failed, using zero fluxes as fallback for all seeds.\n');
        initialSeeds = zeros(length(model.rxns), numSeeds);
        sample.warmUpPoints = initialSeeds;
    end

    % Debugging output to verify warm-up points initialization
    fprintf('Initial warm-up points generated: %d samples with %d reactions each.\n', size(sample.warmUpPoints, 2), size(sample.warmUpPoints, 1));

    % Initialize sample.points with warmUpPoints for ADSB
    omegaSize = rank(sample.warmUpPoints);
    nDim = max([sample.populationScale * omegaSize, 3]); % At least three points required

    % Define the number of points per chain and total chains
    sample.pointsPerChain = nDim;
    sample.numChains = ceil(sample.numSamples / sample.pointsPerChain);

    % Debugging output to verify points per chain and number of chains
    fprintf('Number of chains: %d, Points per chain: %d\n', sample.numChains, sample.pointsPerChain);

    % Initialize points as a 3D array for ADSB
    numReactions = length(model.rxns);
    sample.points = zeros(numReactions, sample.pointsPerChain, sample.numChains);

    % Initialize warm-up points in the 3D structure
    for chainIdx = 1:sample.numChains
        sample.points(:, :, chainIdx) = sample.warmUpPoints(:, 1:sample.pointsPerChain);
    end

    % Run ADSB Sampling using ADSB_v2
    fprintf('Running ADSB_v2 sampling...\n');
    [sample.points, sample.samplingTime] = ADSB_v2(sample);

    % Debugging output to check if sample.points has been populated
    if isempty(sample.points)
        fprintf('Warning: sample.points is empty after ADSB_v2.\n');
    else
        fprintf('Sampling completed with %d samples.\n', size(sample.points, 2) * sample.numChains);
    end

    % Reshape the 3D points matrix to a 2D matrix for final output
    [n, m, p] = size(sample.points);
    sample.points = reshape(sample.points, n, m * p);

    % Verify S*v = 0 for each sample
    tol = 1e-6;  % Tolerance for determining feasibility
    for i = 1:sample.numSamples
        v = sample.points(:, i);  % Extract sample flux vector
        Sv = model.S * v;  % Calculate S*v

        if any(abs(Sv) > tol)
            fprintf('Warning: Sample %d does not satisfy S*v = 0 within tolerance. Max deviation: %e\n', i, max(abs(Sv)));
        end
    end

    fprintf('Sampling completed.\n');
end
