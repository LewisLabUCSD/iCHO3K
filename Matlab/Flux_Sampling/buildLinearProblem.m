function LPproblem = buildLinearProblem(model)
% Builds LP structure problem
%
% Structure has the following form
%
%    min      ~& c^T v \\
%    s.t.     ~& S v = 0 \\
%             ~& lb \leq v \leq ub \\
%
% USAGE:
%              LPproblem = buildLinearProblem(model)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x 1` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%
% OUTPUT:
%              LPproblem (structure):   LP problem structure
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Create new model structure
LPproblem.b = zeros(size(model.S,1),1);
[LPproblem.A,LPproblem.lb,LPproblem.ub] = deal(sparse(model.S),model.lb,model.ub);
LPproblem.rev = 1*((model.lb<0)&(model.ub>0));
LPproblem.c   = zeros(size(LPproblem.ub));

% Constraint sense definition
m = size(model.S,1);
LPproblem.csense = blanks(m);
for ix = 1:m
    LPproblem.csense(ix) = 'E';
end

% Optimization sense (minimization)
LPproblem.osense = 1;
