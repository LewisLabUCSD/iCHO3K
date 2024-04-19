function exRxns = find_ex_rxns(model)

% Note: this function identifies all source or sink reactions, not just those
% exchanging metabolites into and out of the extracellular space

% Find indices of out-only reactions
indOutRxns = sum(model.S > 0) > 0 & sum(model.S < 0) == 0;

% Find indices of in-only reactions
indInRxns = sum(model.S < 0) > 0 & sum(model.S > 0) == 0;

% Define exchange reactions as the union of out- and in-only
indExRxns = indOutRxns | indInRxns;

exRxns = model.rxns(indExRxns);