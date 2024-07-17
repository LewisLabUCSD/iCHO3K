function extracted_models = extract_mCADRE_models(model,ubiData,cell_line,phase, protected_reactions, checkFunctionality)
%%
%   EXTRACT_MCADRE_MODELS extracts context-specific models using the mCADRE
%   algorithm for a given flux-consistent parent genome-scale model MODEL
%   and expression data contained within UBIDATA. 
%
%   INPUTS
%
%   MODEL:                  This is the flux-consistent model in COBRA
%                           format. The model should not contain any
%                           blocked reactions. The lower and upper bounds
%                           of metabolic functionalities must be
%                           pre-specified so that they are quantitatively
%                           protected in the extracted models. (See
%                           Gopalakrishnan et al, 20xx for details)
%
%   UBIDATA:                Structure containing ubiquity scores with the
%                           following fields:
%                           RXNS:       List of reactions. Must be a
%                                       superset of reactions in MODEL.
%                           UBISCORES:  M-by-N matrix of ubiquity scores.
%                                       Elements represent the ubiquity
%                                       score of the mth reaction in nth
%                                       condition.
%
%   CELL_LINE:              Cell line in character string format
%
%   PHASE:                  Process phase in character string format
%
%
%   OUTPUTS
%
%   EXTRACTED_MODELS is a structure with the following fields:
%
%   CELL_LINE:              Cell line in character string format
%
%   PHASE:                  Process phase in character string format
%
%   BASE_MODEL:             The base model from which the context-specific
%                           models are extracted. Same as MODEL.
%
%   RED_MODELS:             Structure with the following fields:
%                           RXNS:           List of reactions in UBIDATA
%                           RETAINED_RXNS:  Matrix of 0s and 1s. Rows
%                                           represent reactions in RXNS,
%                                           columns represent the different
%                                           conditions.

%%  Author: SARATRAM GOPALAKRISHNAN
%   Place:  4A15, BRF-II, UC-San Diego
%   Date: 05/16/2022
%
%% Main Code
extracted_models.cell_line = cell_line;
extracted_models.phase = phase;
extracted_models.base_model = model;

ubiScores = ubiData.ubiScores;

rxns = ubiData.rxns;
ubiScores = ubiScores(ismember(rxns,model.rxns),:);

nm = size(ubiScores,2);
retr = zeros(length(rxns),nm);

environment = getEnvironment();
confScores = ubiData.confidenceScores;
for i = 1:nm
    restoreEnvironment(environment);
    mx = RMF_mCADRE(model,ubiScores(:,i), confScores, protected_reactions, checkFunctionality);
    retr(:,i) = double(ismember(rxns,mx.rxns));
end

extracted_models.red_models.rxns = rxns;
extracted_models.red_models.retained_rxns = retr;

end
