function [PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, salvageCheck, C_H_genes, method)

% Inputs
% - model: original generic model
% - G: list of genes in expression data
% - U: ubiquity scores corresponding to genes in G
% - confidenceScores: literature-based evidence for generic model reactions
% - salvageCheck: option flag for whether to perform functional check for the
%                 nucleotide salvage pathway (1) or not (0)
% - C_H_genes: predefined high confidence reactions (optional)
% - method: 1 = use fastFVA (glpk) to check consistency; 2 = use fastcc & cplex


% Outputs
% - PM: pruned, context-specific model
% - GM: generic model (after removing blocked reactions)
% - C: core reactions in GM
% - NC: non-core reactions in GM
% - Z: reactions with zero expression (i.e., measured zero, not just
%      missing from expression data)
% - model_C: core reactions in the original model (including blocked)
% - pruneTime: total reaction pruning time
% - cRes: result of model checks (consistency/function)
%       - vs. +: reaction r removed from generic model or not
%       1 vs. 2: reaction r had zero or non-zero expression evidence
%       -x.y: removal of reaction r corresponded with removal of y (num.) total
%            core reactions
%       +x.1 vs. x.0: precursor production possible after removal of 
%                      reaction r or not
%       3: removal of reaction r by itself prevented production of required
%          metabolites (therefore was not removed)

if nargin < 7
    method = 1; % fastFVA
end

if nargin < 6
    C_H_genes = [];
end

if nargin < 5
    salvageCheck = 1;
end

%% Generate order for reaction removal

% Gene ubiquity scores are converted to reaction expression evidence to
% define the core (C) and non-core (NC) reaction sets. Inactive reactions
% are identified and removed from the global model to produce the generic
% model (GM) for subsequent pruning. Non-core reactions are ordered first
% by expression and then by connectivity evidence to give the list P. Any
% reactions with zero expression (i.e., associated, but non-expressed
% genes) are also listed in the vector Z.

display('Processing inputs and ranking reactions...')

[GM, C, NC, P, Z, model_C] = ...
    rank_reactions(model, G, U, confidenceScores, C_H_genes, method);

%% Define inputs to the model pruning step

% Define core vs. non-core ratio threshold for removing reactions
eta = 1/3;

% Check functionality of generic model
changeCobraSolver('glpk');
load('precursorMets');

genericStatus = check_model_function(GM, 'requiredMets', precursorMets);

if genericStatus
display('Generic model passed precursor metabolites test');

%% If generic functionality test is passed, prune reactions
    display('Pruning reactions...')
    t0 = clock;
    
    [PM, cRes] = prune_model(GM, P, C, Z, eta, precursorMets, salvageCheck, ...
        method);
    
    pruneTime = etime(clock,t0);
else
    display('Generic model failed precursor metabolites test')
end