function vFeasible = findNearestFeasiblePoint_ADSB(model,vInfeasible,Nint,vTol)
% Estimate closest feasible point from a query point
%
% Structure has the following MIQP form
%
%    min      ~& ||delta||_2 \\
%    s.t.     ~& S vFeasible = 0 \\
%             ~& vFeasible - vInfeasible = delta \\
%             ~& vFeasible is loopless (see buildLooplessProblem.m description)
%             ~& lb \leq vFeasible \leq ub \\
%             ~& 0 \leq delta \\
%
% USAGE:
%              vFeasible = findNearestFeasiblePoint(model, vInfeasible, Nint, vTol)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x 1` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * intRxns - `ni x 1` Matrix with indexes of internal rxns
%                                    * nRxns - Number of reactions in the model (double)
%                                    * warmUpPoints - `n x k` Matrix of warmup points 
%              vInfeasible:  n x 1 infeasible flux vector
%              Nint:  ni x p Null space basis of Sint
%
% OPTIONAL INPUTS:
%              vTol: Numerical flux tolerance  (default 1e-8)
%
% OUTPUT:
%              vFeasible:   n x 1 feasible flux vector
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Handle missing inputs
if (nargin<4); vTol = 1e-8; end

% Assign model parameters
nRxns   = model.numRxns;
intRxns = numel(model.intRxns);
K       = 1e3;                                                             % large positive constant

% Assign previous structure
MIQPmodel = buildLooplessProblem(model,Nint);

% Add slacks for flux vector
MIQPmodel.A = sparse([MIQPmodel.A,zeros(size(MIQPmodel.A,1),nRxns);...
                    eye(nRxns),zeros(nRxns,2*intRxns),-eye(nRxns)]);
MIQPmodel.c = zeros(size(MIQPmodel.A,2),1);
MIQPmodel.F = sparse([zeros(nRxns+2*intRxns,2*nRxns+2*intRxns);zeros(nRxns,nRxns+2*intRxns),eye(nRxns)]);

% Update remaining fields
MIQPmodel.b       = [MIQPmodel.b;vInfeasible];
MIQPmodel.lb      = [MIQPmodel.lb;-K*ones(nRxns,1)];
MIQPmodel.ub      = [MIQPmodel.ub;K*ones(nRxns,1)];
MIQPmodel.vartype = [MIQPmodel.vartype,repmat('C',1,nRxns)];
MIQPmodel.csense  = [MIQPmodel.csense,repmat('E',1,nRxns)];
MIQPmodel.osense  = 1;

% Find closest point to the query
sol = solveCobraMIQP(MIQPmodel);
if strcmp(sol.stat,'OPTIMAL')
    sol.full(abs(sol.full)<vTol) = 0;
    vFeasible = sol.full(1:nRxns);

    % Find random interior point instead
else
    vFeasible = model.warmUpPoints(:,randi(size(model.warmUpPoints,2)));
end
