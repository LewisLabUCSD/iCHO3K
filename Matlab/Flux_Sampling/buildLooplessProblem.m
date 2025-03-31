function looplessProblem = buildLooplessProblem(model,Nint)
% Builds MILP (loopless) structure problem
%
% Structure has the following form
%
%    min      ~& c^T v \\
%    s.t.     ~& S v = 0 \\
%             ~& N^T g = 0 \\
%             ~& -1000*a_i + (1-a_i) \leq g_i \leq -a_i + 1000*(1-a_i) \\
%             ~& lb*(1-a_i) \leq v_i \leq ub*a_i \\
%             ~& a_i binary
%
% USAGE:
%              looplessProblem = buildLooplessProblem(model, Nint)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x n` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * intRxns - `ni x 1` Matrix with indexes of internal rxns
%                                    * exchRxns - `ne x 1` Matrix with indexes of exchange rxns
%              Nint:  ni x p Null space basis of Sint
%
% OUTPUT:
%              looplessProblem (structure):   MILP problem structure
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Determine number of internal reactions and define bounds
m      = numel(model.intRxns);
[l,n]  = size(model.S);
lb_int = model.lb(model.intRxns);
ub_int = model.ub(model.intRxns);

% Determine constrained internal reactions
M = eye(n);
M(model.exchRxns,:) = [];

% Define model objective and osense
looplessProblem.c      = zeros(n+2*m,1);
looplessProblem.osense = 1;

% LHS matrix definition
K = 1e3;
p = size(Nint,1);
looplessProblem.A = sparse([model.S,zeros(l,2*m);...
                        zeros(p,n),Nint,zeros(p,m);...
                        zeros(m,n),eye(m),(1+K)*eye(m);...
                        zeros(m,n),-eye(m),-(1+K)*eye(m);...
                        M,zeros(m),-diag(ub_int);...
                        -M,zeros(m),-diag(lb_int)]);
                   
% RHS definition
looplessProblem.b = [zeros(l+p,1);K*ones(m,1);-ones(m,1);zeros(m,1);-lb_int];

% Constraint sense definition
looplessProblem.csense = blanks(l+p+4*m);
for ix = 1:l+p+4*m
    if ix <= l+p
        looplessProblem.csense(ix) = 'E';
    else
        looplessProblem.csense(ix) = 'L';
    end
end

% Assignation of variable types
looplessProblem.vartype = blanks(n+2*m);
for ix = 1:n+2*m
    if ix <= n+m
        looplessProblem.vartype(ix) = 'C';
    else
        looplessProblem.vartype(ix) = 'B';
    end
end

% Bounds, reversibilities and initial solution field definition 
looplessProblem.lb  = [model.lb;-K*ones(m,1);zeros(m,1)];
looplessProblem.ub  = [model.ub;K*ones(m,1);ones(m,1)];
looplessProblem.rev = 1*((model.lb<0)&(model.ub>0));
looplessProblem.x0  = [];
