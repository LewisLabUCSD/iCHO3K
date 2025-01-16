function [model,rxnList] = parseInternalRxns(model,rxnList)
% Parse reactions from a model
%
% Reorganizes rxns from a model into internal and external
%
% USAGE:
%              model = parseInternalRxns(model)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S - `m x n` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
% OPTIONAL INPUTS:
%              rxnList:  k x 1 Reaction list (cell array)
%
% OUTPUT:
%              model (structure):  COBRA model with rxns reorganized
%
% OPTIONAL OUTPUT:
%              rxnList:  k x 1 Updated rxnList
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

if (nargin<2); rxnList = []; end

% Find internal and exchange reactions and re-arrange stoichiometric matrix
intRxns  = find(sum(full(model.S~=0))>1);
exchRxns = find(sum(full(model.S~=0))<=1);
model.S  = model.S(:,[intRxns,exchRxns]);

% Reorganize reaction fields
model.lb = model.lb([intRxns,exchRxns]);
model.ub = model.ub([intRxns,exchRxns]);
if ~isempty(rxnList)
    rxnList = [rxnList(intRxns);rxnList(exchRxns)];
end
model.intRxns  = 1:numel(intRxns);
exchRxns       = 1:numel(exchRxns);
model.exchRxns = exchRxns + numel(intRxns);
