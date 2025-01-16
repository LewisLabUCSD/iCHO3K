function [model,zeroRxns] = removeBlockedSets(model)
% Remove blocked sets
%
% Removes blocked rxns and mets sets from a model
%
% USAGE:
%              model = removeBlockedSets(model)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S/A - `m x n` Stoichiometric matrix
%                                    * c  - `n x 1` Objective vector
%                                    * b  - `m x 1` RHS vector
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%
% OUTPUT:
%              model (structure):  COBRA model without blocked sets
%
% OPTIONAL OUTPUT:
%              zeroRxns:  k x 1 Index vector of blocked rxns
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

zeroRxns = find((model.ub==0)&(model.lb==0));
if ~isempty(zeroRxns)
    if isfield(model,'A')
        model.A(:,zeroRxns) = [];
        model.c             = zeros(size(model.A,2),1);
        
        % Find orphan metabolites and remove
        model.A(all(model.A==0,2),:) = [];
        
        % Re-define other quantities
        model.b = zeros(size(model.A,1),1);
    end
    
    if isfield(model,'S')
        model.S(:,zeroRxns) = [];
        model.c = zeros(size(model.S,2),1);
        
        % Find orphan metabolites and remove
        model.S(all(model.S==0,2),:) = [];
        
        % Re-define other quantities
        model.b = zeros(size(model.S,1),1);
    end
    
    % Remove rxn-related fields
    model.lb(zeroRxns) = [];
    model.ub(zeroRxns) = [];
end
model.rev = (model.lb<0)&(model.ub>0);
