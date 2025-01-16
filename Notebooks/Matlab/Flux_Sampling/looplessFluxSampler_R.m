function sample = looplessFluxSampler_R(model,options)
% Loopless Flux Sampler
%
% Performs random sampling of the loopless space of metabolic models
%
% USAGE:
%              sample = looplessFluxSampler(model, options)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x 1` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * rxns - `n x 1` rxn identifiers (cell array)
%              options (structure):  (the following fields are required, others can be supplied - see below)
%                                    * numSamples - number of points (double)
%
% OPTIONAL INPUTS:
%              options (structure):  (the following fields are optional)
%                                    * numDiscarded - Burn-in (double) (only used in ll-ACHRB) (default 0)
%                                    * stepsPerPoint - Thinning or number of steps per effective point (double). In ADSB, it refers to the expected number of times a point is moved
%                                    * algorithm - 'ADSB' (default) or 'll_ACHRB' or 'HR'
%                                    * loopless - Loop removal flag, true (default) or false
%                                    * vTol - Numerical flux tolerance (default 1e-8)
%                                    * parallelFlag - Parallel sampling option, false (default) or true
%                                    * warmUpFlag - Warmup option, true (default) or false
%                                    * numCores - Number of cores for parallel sampling (double), empty (default)
%                                    * diagnostics - MCMC diagnostics flag, true (default) or false
%                                    * populationScale - Size of current set relative to the feasible space dimension (integer => 1), 3 (default)
%                                    * points - Points from a previous sampling, empty (default)
%
% OUTPUT:
%              sample (structure):   sampling structure containing #numSamples random points
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Check inputs
if (nargin<2)
    fprintf('Not enough input arguments.\n');
    return;
else
    % Check if COBRA toolbox is in the path
    if ~exist('initCobraToolbox','file')
        fprintf('COBRA Toolbox is not in the path. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
        return;
    end
    
    % Number of samples
    if isfield(options,'numSamples'); sample.numSamples = options.numSamples;
    else fprintf('Field numSamples has not been defined.'); return; end
    
    % Samples discarded or burn-in (default 0)
    if isfield(options,'numDiscarded'); sample.numDiscarded = options.numDiscarded;
    else sample.numDiscarded = 0; end
    
    % Steps per point (thinning in ACHRB) or number of expected moves (ADSB) (default 1e2)
    if isfield(options,'stepsPerPoint'); sample.stepsPerPoint = options.stepsPerPoint;
    else sample.stepsPerPoint = 1e2; end
    
    % Sampling algorithm 'ADSB' (default) or 'll_ACHRB' or 'HR'
    if isfield(options,'algorithm'); sample.algorithm = options.algorithm;
    else sample.algorithm = 'ADSB'; end
    
    % Sampling mode (default loopless true or 1)
    if isfield(options,'loopless'); sample.loopless = options.loopless;
    else sample.loopless = 1; end
    
    % Numerical flux tolerance  (default 1e-8)
    if isfield(options,'vTol'); sample.vTol = options.vTol;
    else sample.vTol = 1e-8; end
    
    % Parallel sampling option (default false or 0)
    if isfield(options,'parallelFlag'); sample.parallelFlag = options.parallelFlag;
    else sample.parallelFlag = 0; end
    
    % Number of cores for parallel sampling (default empty) (at least two cores for parallel computations
    if isfield(options,'numCores'); sample.numCores = max([2,options.numCores]);
    else sample.numCores = []; end
    
    % Warmup flag (default true or 1)
    if isfield(options,'warmUpFlag'); sample.warmUpFlag = options.warmUpFlag;
    else sample.warmUpFlag = 1; end
    
    % Run MCMC diagnostics
    if isfield(options,'diagnostics'); sample.diagnostics = options.diagnostics;
    else sample.diagnostics = 1; end
    
    % Fold number of particles in the population relative to dim(Omega) (only for ADSB, default 1)
    if strcmp(sample.algorithm,'ADSB') && isfield(options,'populationScale'); sample.populationScale = options.populationScale;
    else sample.populationScale = 3; end
    
    % Restart sampler from previous results
    if ~isfield(model,'points'); sample.points = [];
    else sample.points = model.points; end
end

% Define solver parameters
changeCobraSolverParams('MILP','intTol',1e-9);
changeCobraSolverParams('MILP','relMipGapTol',1e-12);
changeCobraSolverParams('LP','optTol',1e-9);
changeCobraSolverParams('LP','feasTol',1e-9);

% Initialize temporal model structure
tempModel = model;
[tempModel,rxnList] = parseInternalRxns(tempModel,model.rxns);
fprintf('Model loaded: initial model contains %d rxns and %d mets.\n',numel(model.rxns),numel(model.mets));

%% I. Model reduction
% Time model pre-processing and sampling preparation
t0 = cputime;

% Reduce model by removing singletons. First, build LP structure
LPproblem = buildLinearProblem(tempModel);

% Run FVA to remove singletons
[tempModel.lb,tempModel.ub] = generalFVA(LPproblem,sample.vTol);
[tempModel,zeroRxns] = removeBlockedSets(tempModel);
rxnList(zeroRxns)    = [];
[tempModel,rxnList]  = parseInternalRxns(tempModel,rxnList);

% Run ll-FVA to remove rxns involved in infeasible loops active
if (sample.loopless)
    
    % Get a sparse null-space matrix of internal rxns
    Nint = fast_snp(tempModel.S(:,tempModel.intRxns),-1e2*(tempModel.lb(tempModel.intRxns)<0),...
        1e2*(tempModel.ub(tempModel.intRxns)>0),sample.vTol);
    if ~isempty(Nint)
        looplessProblem = buildLooplessProblem(tempModel,normMatrixEntries(Nint)');
        [tempModel.lb,tempModel.ub] = generalFVA(looplessProblem,sample.vTol,'loopless');
        
        % If there are no possible active loops, then turn-off loopless option and
        % run conventional FVA (not ll-FVA)
    else
        sample.loopless = 0;
        LPproblem       = buildLinearProblem(tempModel);
        [tempModel.lb,tempModel.ub] = generalFVA(LPproblem,sample.vTol);
    end
    [tempModel,zeroRxns] = removeBlockedSets(tempModel);
    rxnList(zeroRxns)    = [];
    [tempModel,rxnList]  = parseInternalRxns(tempModel,rxnList);
    
    % Get an updated sparse null-space matrix of internal rxns (if relevant)
    if sample.loopless
        Nint = fast_snp(tempModel.S(:,tempModel.intRxns),-1e2*(tempModel.lb(tempModel.intRxns)<0),...
            1e2*(tempModel.ub(tempModel.intRxns)>0),sample.vTol);
        if ~isempty(Nint); Nint = normMatrixEntries(Nint);
        else sample.loopless = 0; end
    end
    if sample.loopless; fprintf('The model has %d infeasible loop law(s) potentially active.\n',size(Nint,2));
    else fprintf('The model does not contain infeasible loops active.\n'); end
end
fprintf('Model reduced: pre-processed model has %d rxns.\n',numel(rxnList));

%% II. Preparation phase
% Define fields related to the model in sample structure
sample.rxns     = rxnList;
sample.numRxns  = numel(tempModel.lb);
sample.S        = tempModel.S;
sample.lb       = tempModel.lb;
sample.ub       = tempModel.ub;
sample.intRxns  = tempModel.intRxns;
sample.exchRxns = tempModel.exchRxns;
sample.uTol     = 1e-10;                                                   % Direction tolerance
sample.bTol     = 1e-10;                                                   % Minimum allowed distance to closest constraint
if sample.loopless
    sample.Nint       = Nint;
    sample.nnzEntries = sum(Nint~=0);
    sample.loopMatrix = sparse(sign(Nint));
    
    % Define loopless feasibility oracle
    sample.isFeasible = @(points) looplessCheck(points,sample);
end

% Define fxn to ensure samples are kept within the bounds
sample.keepWithinBounds = @(points) bringToBoundary(points,sample.lb,sample.ub);

% Generate initial (loopless) flux seeds
sample.warmUpPoints = generateLooplessFluxSeeds(sample);

% Calculate centroid
centroid = mean(sample.warmUpPoints,2);
centroid(abs(centroid)<sample.vTol) = 0;

% Check feasibility of the initial point. If centroid infeasible, find closest feasible point
if sample.loopless && ~sample.isFeasible(centroid)
    sample.centroid = findNearestFeasiblePoint(sample,centroid,Nint',sample.vTol);
else
    sample.centroid = centroid;
end

% Generate seeds for running ADSB and building current set
if strcmp(sample.algorithm,'ADSB')
    
    % Determine size of the feasible space
    omegaSize = rank(sample.warmUpPoints);
    if (sample.populationScale == 1); nDim = sample.populationScale*omegaSize+1;
    else nDim = sample.populationScale*omegaSize; end
    
    % Define number of chains and points per chain (at least three points are required)
    minPoints             = 3;
    sample.pointsPerChain = max([nDim,minPoints]);
    sample.numChains      = ceil(sample.numSamples/sample.pointsPerChain);
    
    % Generate seed points for the current set: Run ll-ACHRB
    if isempty(sample.points)
        nWarmups   = 2e4;
        maxPoints  = 2e6;
        cWeight    = size(sample.warmUpPoints,2);
        cCentroid  = [];
        prevPoints = [];
        while (nWarmups < maxPoints)
            [vPoints,~,sample.centroid,cWeight] = ll_ACHRB(sample,nWarmups,0,1,1,cWeight);
            if isempty(vPoints); return; end
            
            % Calculate centroid trace evolution
            prevPoints = [prevPoints,vPoints];
            
            % Check if warmup flag is on
            if ~sample.warmUpFlag; break; end
            
            % Calculate current centroid
            cCentroid  = [cCentroid,bsxfun(@rdivide,cumsum(vPoints,2),size(cCentroid,2)+1:size(cCentroid,2)+size(vPoints,2))];
            Rfactor    = psrf(cCentroid');
            
            % Check the convergence of the reduction factor
            if (abs(median(Rfactor(isfinite(Rfactor)))-1) <= .1); break;
            else nWarmups = 2*nWarmups; end
        end
    else
        % Use previous points in the current set
        prevPoints = sample.points;
    end
    
    % Define current sets
    if isempty(sample.points)
        sample.points = zeros(sample.numRxns,nDim,sample.numChains);
        for ix = 1:sample.numChains
            sample.points(:,:,ix) = prevPoints(:,sampleSet(size(prevPoints,2),nDim,'false'));
        end
    else
        extraPoints = sample.pointsPerChain*sample.numChains-sample.numSamples;
        if (extraPoints ~= 0)
            extraPoints = prevPoints(:,sampleSet(size(prevPoints,2),extraPoints,'false'));
            prevPoints   = [prevPoints,extraPoints];
        end
        prevPoints = reshape(prevPoints,sample.numRxns,nDim,sample.numChains);
        sample.points = prevPoints;
    end
end

% Report ADSB initial stats
if strcmp(sample.algorithm,'ADSB')
    fprintf('%d points spanning the feasible space have been incorporated to the current set.\n',sample.pointsPerChain);
else
    fprintf('%d warmup points spanning the feasible space have been defined.\n',size(sample.warmUpPoints,2));
end
sample.prepTime = (cputime-t0)/60;
clearvars -except sample

%% III. Sampling
% Run appropriate sampler in either single core or parallel mode
% Single core option
if ~sample.parallelFlag
    fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
    if strcmp(sample.algorithm,'ADSB')
        [sample.points,sample.samplingTime] = ADSB(sample);
        [n,m,p] = size(sample.points);
        sample.points = reshape(sample.points,n,m*p);
        sample.points = sample.points(:,1:sample.numSamples);

    elseif strcmp(sample.algorithm,'ll_ACHRB')
        [sample.points,sample.samplingTime,sample.centroid] = ll_ACHRB(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint);
        
    elseif strcmp(sample.algorithm,'HR')
        [sample.points,sample.samplingTime,sample.rejectionRate] = HitAndRun(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint);
    end

    if isempty(sample.points); return; end
    
    % Multiple core option
else
    fprintf('Sampling in progress (multiple cores)...\n');
        
    % Delete active workers
    if ~isempty(gcp('nocreate')); delete(gcp('nocreate')); end        
   
    % Try to use the number of workers requested. If not posible, use the number of workers allocated to the 'local' profile
    if ~isempty(sample.numCores)
        try
            parpool(sample.numCores);
        catch
            parWorkers = parpool('local');
            sample.numCores = parWorkers.NumWorkers;
        end
    else
        parWorkers = parpool('local');
        sample.numCores = parWorkers.NumWorkers;
    end
    
    % Choose appropriate sampling algorithm
    if strcmp(sample.algorithm,'ADSB')
        
        % Replicate sample structure for parallel sampling
        chainsPerWorker = fix(sample.numChains/sample.numCores);
        chainsPerWorker = chainsPerWorker*ones(1,sample.numCores-1);
        chainsPerWorker = [chainsPerWorker,sample.numChains-sum(chainsPerWorker)];
        
        % Assign appropriate number of structures to each core
        idx1 = 0; idx2 = 0;
        for ix = 1:sample.numCores            
            idx1 = idx1+1;                                 % Update indices
            idx2 = idx2+chainsPerWorker(ix);
            samples{ix}           = sample;
            samples{ix}.numChains = chainsPerWorker(ix);
            samples{ix}.points    = sample.points(:,:,idx1:idx2);
            idx1 = idx2;                           % Update starting points
        end

        % Perform parallel sampling using ADSB
        parfor workerIdx = 1:sample.numCores
            rng('shuffle');
            [samples{workerIdx}.points,samples{workerIdx}.samplingTime] = ADSB(samples{workerIdx},workerIdx);
        end
        delete(gcp('nocreate'));

        % Build definitive structure
        sample.points       = [];
        sample.samplingTime = 0;
        for ix = 1:sample.numCores
            sample.points = [sample.points,reshape(samples{ix}.points,...
                                    sample.numRxns,sample.pointsPerChain*chainsPerWorker(ix))];
            sample.samplingTime = max([sample.samplingTime,samples{ix}.samplingTime]);
        end

        % Fill with the exact number of samples
        sample.points = sample.points(:,1:sample.numSamples);

    elseif strcmp(sample.algorithm,'ll_ACHRB') || strcmp(sample.algorithm,'HR')
        
        % Perform parallel sampling using ll-ACHRB or HR
        numPointsPerCore              = round(sample.numSamples/sample.numCores);
        samplingTime{sample.numCores} = 0;        
        points{sample.numCores}       = [];
        if strcmp(sample.algorithm,'ll_ACHRB')
            centroid{sample.numCores} = [];
            parfor workerIdx = 1:sample.numCores
                rng('shuffle');
                [points{workerIdx},samplingTime{workerIdx},centroid{workerIdx}] = ll_ACHRB(sample,numPointsPerCore,sample.numDiscarded,sample.stepsPerPoint,1,workerIdx);
            end

        elseif strcmp(sample.algorithm,'HR')   
            parfor workerIdx = 1:sample.numCores
                rng('shuffle');
                [points{workerIdx},samplingTime{workerIdx}] = HitAndRun(sample,numPointsPerCore,sample.numDiscarded,sample.stepsPerPoint,1,workerIdx);
            end
        end
        delete(gcp('nocreate'));

        % Build final sample structure
        sample.samplingTime = 0;
        sample.points       = [];
        sample.centroid     = [];
        for ixCore = 1:sample.numCores
            sample.points       = [sample.points,points{ixCore}];
            sample.samplingTime = max([sample.samplingTime,samplingTime{ixCore}]);
            if strcmp(sample.algorithm,'ll_ACHRB')
                sample.centroid = [sample.centroid,centroid{ixCore}];
            end
        end

        % Fill with the exact number of samples
        sample.points = sample.points(:,1:sample.numSamples);
        if strcmp(sample.algorithm,'HR')
            sample.centroid = sum(sample.points,2);
        else
            sample.centroid = sum(sample.centroid,2);
        end
    end
end
fprintf('---------------------------\nSampling finalized.\n');

% Clear variables and remove auxiliary fields
clearvars -except sample
sample = rmfield(sample,'warmUpPoints');

%%  IV. Run sampling and MCMC diagnostics
if sample.diagnostics
    fprintf('Starting MCMC diagnostics...\n');
    
    % Split samples resulting from the chain in 10 segments for
    % convergence analysis
    numChainsForDiagnostics      = 10;
    pointsPerChainForDiagnostics = fix(sample.numSamples/numChainsForDiagnostics);
    sample.points = reshape(sample.points,sample.numRxns,pointsPerChainForDiagnostics,numChainsForDiagnostics);
    
    % Calculate potential scale reduction statistics and recover original chain
    sample.points = permute(sample.points,[2,1,3]);
    [sample.R,sample.Rint,sample.Neff,sample.tau,sample.thin] = psrf(sample.points);
    sample.points = permute(sample.points,[2,3,1]);
    sample.points = reshape(sample.points,sample.numRxns,pointsPerChainForDiagnostics*numChainsForDiagnostics);
    sample.points = sample.points(:,1:sample.numSamples);
    
    % Calculate chains statistics
    sample.mu     = mean(sample.points,2);
    sample.sigma  = std(sample.points')';
    sample.median = median(sample.points,2);
    sample.IQR    = iqr(sample.points,2);
    
    % Check if there are exclusive topological modes
    if (sample.loopless)
        
        % Determine unique flux patterns and continue sampling (if necessary)
        fluxPattern = (~sample.sigma)&(~sample.mu);
        if any(fluxPattern)
            uniquePatterns = unique(fluxPattern','rows');
            fprintf('%d blocked reaction(s) detected in the final sample.\n',size(uniquePatterns,1));
        end
    end
end