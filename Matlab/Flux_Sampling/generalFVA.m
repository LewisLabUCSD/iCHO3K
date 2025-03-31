function [minFlux,maxFlux,Vmin,Vmax] = generalFVA(optProblem,vTol,type)
% Performs General Flux Variability Analysis (allowing or blocking
% infeasible loops)
%
% Structure has the following MILP form
%
%    min/max  ~& v_i \\
%    s.t.     ~& S v = 0 \\
%             ~& lb \leq v \leq ub \\
%             ~& v is loopless if type is 'loopless'
%              (see buildLooplessProblem.m for loopless MILP formulation description)
%
% USAGE:
%              [minFlux,maxFlux] = generalFVA(optProblem)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x 1` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * intRxns - `ni x 1` Matrix with indexes of internal rxns
%                                    * nRxns - Number of reactions in the model (double)
%                                    * rev - Boolean vector indicating reversible reactions
%              Nint:  ni x p Null space basis of Sint
%
% OPTIONAL INPUTS:
%              vTol: Numerical flux tolerance  (default 1e-8)
%              type: FVA mode, 'loopless' (default) or 'conventional'
%
% OUTPUT:
%              minFlux:   n x 1 flux vector of minimum fluxes for each rxn
%              maxFlux:   n x 1 flux vector of maximum fluxes for each rxn
%
% OPTIONAL OUTPUTS:
%              Vmin:   n x n Flux matrix with the solutions associated to minFlux
%              Vmax:   n x n Flux matrix with the solutions associated to maxFlux
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Define parameters and type of FVA optimization
if (nargin<2)
    vTol       = 1e-8;
    allowLoops = 1;
elseif (nargin<3)
    allowLoops = 1;
else
    if strcmp(type,'loopless'); allowLoops = 0;
    else allowLoops = 1; end
end

% Pre-allocate memory for the subsequent optimizations
Vmin = zeros(size(optProblem.rev,1));
Vmax = zeros(size(optProblem.rev,1));

for ix = 1:numel(optProblem.rev)
    
    optProblem.c(ix) = 1;
    
    % Minimization problem
    optProblem.osense = 1;
    if allowLoops; solMin = solveCobraLP(optProblem);
    else solMin = solveCobraMILP(optProblem); end
    
    % Maximization problem
    optProblem.osense = -1;
    if allowLoops; solMax = solveCobraLP(optProblem);
    else solMax = solveCobraMILP(optProblem); end
    
    % Assign optimal values
    Vmin(:,ix) = solMin.full(1:numel(optProblem.rev));
    Vmax(:,ix) = solMax.full(1:numel(optProblem.rev));
    
    % Reset objective
    optProblem.c(ix) = 0;
end

% Solutions from the pre-processing step have fluxes greater tolerance
Vmin(abs(Vmin)<vTol) = 0;
Vmax(abs(Vmax)<vTol) = 0;

% Assign minima and maxima
minFlux = diag(Vmin);
maxFlux = diag(Vmax);

% Ensure that the solutions found are contained in the feasible region
Vmax = bringToBoundary(Vmax,minFlux,maxFlux);
Vmin = bringToBoundary(Vmin,minFlux,maxFlux);
