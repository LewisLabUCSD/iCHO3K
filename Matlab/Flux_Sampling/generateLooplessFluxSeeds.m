function warmupPoints = generateLooplessFluxSeeds(sample,vTol)
% Generates a set of (loopless) flux points
%
% Uses the warmup method of ll-ACHRB to generate a (loopless) sample of
% flux solutions
%
% USAGE:
%              warmupPoints = warmupLooplessACHRB(sample)
%
% INPUTS:
%              sample (structure):    (the following fields are required - others can be supplied)
%                                    * S - `m x n` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * numRxns - Number of reactions in the model (double)
%                                    * Nint - ni x p Null space basis of Sint
%                                    * intRxns - `ni x 1` Matrix with indexes of internal rxns
%                                    * exchRxns - `ne x 1` Matrix with indexes of exchange rxns
%                                    * loopless - Loop removal option, true (default) or false
%
% OPTIONAL INPUTS:
%              vTol: Numerical flux tolerance  (default 1e-8)
%
% OUTPUT:
%              warmupPoints:  Set of (loopless) flux solutions
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Define number of warmup points
if (nargin<2); vTol = 1e-8; end

% Initialize main variables
warmupPoints = zeros(sample.numRxns,2*sample.numRxns);
LPproblem    = buildLinearProblem(sample);

% Perform loopless warmup
if sample.loopless; looplessProblem = buildLooplessProblem(sample,sample.Nint'); end

% Define parameters for the warmup process
lbRef = LPproblem.lb;
ubRef = LPproblem.ub;
alpha = 0.95;

% Main loop
for counter = 1:sample.numRxns
    
    % Change objective
    LPproblem.c(counter)  = 1;
    
    % a) Solve maximization problem (try first solving the LP problem)
    LPproblem.osense      = -1;    
    LPproblem.ub(counter) = lbRef(counter) + alpha*(ubRef(counter)-lbRef(counter));
    solMax                = solveCobraLP(LPproblem);
    solMax.full(abs(solMax.full)<vTol) = 0;
    LPproblem.ub(counter) = ubRef(counter);
    
    % Check feasibility
    if ~sample.loopless || (sample.loopless && sample.isFeasible(solMax.full))
        warmupPoints(:,2*counter-1) = solMax.full;        
        
    % Run MILP if solution is infeasible
    else
        looplessProblem.c(counter)  = 1;
        looplessProblem.osense      = -1;        
        looplessProblem.ub(counter) = lbRef(counter) + alpha*(ubRef(counter)-lbRef(counter));
        solMax                      = solveCobraMILP(looplessProblem);
        solMax.full(abs(solMax.full)<vTol) = 0;
        warmupPoints(:,2*counter-1) = solMax.full(1:sample.numRxns);
        looplessProblem.ub(counter) = ubRef(counter);
        looplessProblem.c(counter)  = 0;
    end   
        
    % Solve minimization problem (try first solving the LP problem)
    LPproblem.osense      = 1;
    LPproblem.lb(counter) = lbRef(counter) + (1-alpha)*(ubRef(counter)-lbRef(counter));
    solMin                = solveCobraLP(LPproblem);
    solMin.full(abs(solMin.full)<vTol) = 0;
    LPproblem.lb(counter) = lbRef(counter);
    
    % Check feasibility of the LP solution
    if ~sample.loopless || (sample.loopless && sample.isFeasible(solMin.full))
        warmupPoints(:,2*counter) = solMin.full;
        
        % Run MILP problem if solution is infeasible
    else
        looplessProblem.c(counter)      = 1;
        looplessProblem.osense(counter) = 1;
        looplessProblem.lb(counter)     = lbRef(counter) + (1-alpha)*(ubRef(counter)-lbRef(counter));
        solMin                          = solveCobraMILP(looplessProblem);
        solMin.full(abs(solMin.full)<vTol) = 0;
        warmupPoints(:,2*counter-1)     = solMin.full(1:sample.numRxns);
        looplessProblem.lb(counter)     = lbRef(counter);
        looplessProblem.c(counter)      = 0;
    end
    
    % Reset original obj. values
    LPproblem.c(counter) = 0;
end
