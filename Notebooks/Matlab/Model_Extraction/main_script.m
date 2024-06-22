% Step 1: Initialize the COBRA toolbox
initCobraToolbox(0);

path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';
addpath('/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells/Notebooks/Matlab/Standep');

% Step 2: Ensure that the parent genome-scale model is loaded into the
% workspace and assigned the variable "model"
model = genome_scale_model; % Replace"genome_scale_model" with the name of your actual model loaded.
% Note: The genome-scale model must be flux-consistent (having no reactions with zero lower and upper bounds) with the
% extracellular flux measurements already incorporated. If the model has
% not been pre-processed, use the function FluxVariability

% Step 3: Ensure that the ubiquity scores are loaded into the MATLAB
% workspace
UbiData = ubiquity_Scores_variable; %   Replace "ubiquity_scores_variable" with the actual variable from the workspace.

% Step 4: Model extraction using mCADRE
cell_line = 'B';    % Replace with actual cell line name here
phase = 'P';        % Replace with actual phase label here.
i = 1;              % Replace with index of desired replicate. This depends on the number of columns in UbiData.ubiScores
UbiData.ubiScores = UbiData.ubiScores(:,i);
changeCobraSolverParams('LP','feasTol',1e-9);
extracted_models = extract_mCADRE_models(model,UbiData,cell_line,phase);

% Step 5: Retrieve extracted model
r1 = extracted_models.red_models.rxns(any(extracted_models.red_models.retained_rxns,2));
reduced_model = removeRxns(model,model.rxns(~ismember(model.rxns,r1)));
reduced_model = removeUnusedGenes(reduced_model);
