function [inactiveRxns, time, result] = check_model_consistency(model, method, r, deCheck, C)


% This function is designed to quickly identify dead-end reactions in a
% stoichiometric model. The algorithm is largely based on the heuristic
% speed-up to Flux Variability Analysis (FVA) proposed by Jerby et al. [1],
% with modifications to further reduce computation time in Matlab. The
% function can operate independently to report the inactive reactions for
% an entire model, or within a pruning algorithm (e.g., MBA) to examine the
% effect of removing reactions.

% Inputs:
% - model: COBRA model structure
% (optional inputs:)
% - method: parameter specifying whether to use fastFVA (1) or fastcc (2)
% - r: name of reaction to be removed (for model pruning in mCADRE or MBA)
% - deCheck: check for core reactions containing dead end metabolites (only for
%            use with model pruning in MBA)
% - C: list of core reaction names (only for model pruning in MBA)

% Outputs:
% - inactiveRxns: list of IDs corresponding to reactions with 0 mininum and
%                 0 maximum flux
% - time: CPU time required to complete function
% - result: summary indicator of dead-end effects on inactive reactions
%       1: removal of r did not create metabolite dead ends leading to
%          inactivation of core reactions
%       2: removal of r created metabolite dead ends leading to
%          inactivation of core reactions

if nargin < 4
    deCheck = 0;
    C = {};
end

if nargin < 3
    r = [];
end

if nargin < 2
    method = 1;
end


if numel(r)
   % Remove reaction r from the model
    model = removeRxns(model, r);
end
model.c(logical(model.c)) = 0;

inactiveRxns = r;

t0 = clock;
result = 1; % Until proven otherwise, assume that removal of r does not
            % create any metabolite dead ends

%% First check whether any core reactions are blocked by the removal of r.

% Checking for metabolite dead ends is accomplished entirely by matrix
% operations, and is therefore very fast in Matlab. If any core reaction
% contains a dead-end metabolite, the reaction itself will be a dead end.
% This check potentially avoids sequential optimizations, as the function
% can exit if any blocked core reactions are detected.
if deCheck
    deadEnd_C = check_core_deadends(model, C);
else
    deadEnd_C = [];
end

% If no core reactions were found to be blocked based on metabolite dead
% ends, maximize and minimize reactions to identify those with zero flux
% capacity.
if numel(deadEnd_C)
    % Setting inactiveRxns to include dead-end containing reactions will
    % effectively cause the function to exit without checking non-core
    % reactions below; thus, the full list of inactive reactions will not
    % be enumerated
    inactiveRxns = union(deadEnd_C, inactiveRxns);

    % This updates the indicator to report that dead-end-containing
    % reactions were found in the core
    result = 2;

% If the option is specified, fastFVA is used to quickly scan through all
% reactions. **note: may want to include option to use fastFVA with GLPK
else
    inactiveRxns = union(inactiveRxns, find_inactive_rxns(model, method));
end

time = etime(clock,t0);
display(['check_model_consistency time: ',num2str(time, '%1.2f'), ' s'])
