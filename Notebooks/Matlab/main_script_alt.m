% Step 1: Initialize the COBRA toolbox
initCobraToolbox(0);

path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';
addpath('/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells/Notebooks/Matlab/Model_Extraction');

% Step 2: Ensure that the parent genome-scale model is loaded into the
% workspace and assigned the variable "model"
model = readCbModel(fullfile(path, 'Notebooks', 'iCHO3595_unblocked.mat'));
% Note: The genome-scale model must be flux-consistent (having no reactions with zero lower and upper bounds) with the
% extracellular flux measurements already incorporated. If the model has
% not been pre-processed, use the function FluxVariability

% Step 3: Ensure that the ubiquity scores are loaded into the MATLAB
% workspace
matFilePath = fullfile(path, 'Data/Context_specific_models', 'ubiData.mat');
loadedData = load(matFilePath);

% Check if the expected variable exists in the loaded data
if isfield(loadedData, 'ubiData')
    UbiData = loadedData.ubiData;
else
    error('The variable "ubiData" does not exist in the loaded .mat file.');
end

% Verify that UbiData is a structure and contains the field 'ubiScores'
if ~isstruct(UbiData)
    error('UbiData is not a structure.');
end

if ~isfield(UbiData, 'ubiScores')
    error('UbiData does not contain the field "ubiScores".');
end

% Display contents of UbiData to verify
disp('Contents of UbiData:');
disp(UbiData);

% Step 4: Model extraction using mCADRE
sampleConditions = {
    'S1', 'WT_P2'; 'S7', 'WT_P2'; 'S13', 'WT_P2';
    'S2', 'WT_P4'; 'S8', 'WT_P4'; 'S14', 'WT_P4';
    'S3', 'WT_P6'; 'S9', 'WT_P6'; 'S15', 'WT_P6';
    'S4', 'WT_P8'; 'S10', 'WT_P8'; 'S16', 'WT_P8';
    'S5', 'WT_P12'; 'S11', 'WT_P12';
    'S6', 'WT_P14'; 'S12', 'WT_P14'; 'S17', 'WT_P14';
    'S18', 'ZeLa_P4'; 'S23', 'ZeLa_P4'; 'S28', 'ZeLa_P4'; 'S34', 'ZeLa_P4'; 'S39', 'ZeLa_P4';
    'S19', 'ZeLa_P6'; 'S29', 'ZeLa_P6'; 'S35', 'ZeLa_P6'; 'S40', 'ZeLa_P6';
    'S20', 'ZeLa_P8'; 'S24', 'ZeLa_P8'; 'S30', 'ZeLa_P8'; 'S36', 'ZeLa_P8'; 'S41', 'ZeLa_P8';
    'S25', 'ZeLa_P12'; 'S31', 'ZeLa_P12'; 'S42', 'ZeLa_P12';
    'S21', 'ZeLa_P14'; 'S26', 'ZeLa_P14'; 'S32', 'ZeLa_P14'; 'S37', 'ZeLa_P14'; 'S43', 'ZeLa_P14';
    'S22', 'ZeLa_P2'; 'S27', 'ZeLa_P2'; 'S33', 'ZeLa_P2'; 'S38', 'ZeLa_P2';
};

for i = 1:length(UbiData.ubiScores(1,:))
    % Determine the corresponding cell line and phase
    sample = sampleConditions{i, 1};
    condition = sampleConditions{i, 2};
    splitCondition = strsplit(condition, '_');
    cell_line = splitCondition{1};
    phase = splitCondition{2};

    % Update the ubiquity scores for the current sample
    currentUbiScores = UbiData;
    currentUbiScores.ubiScores = UbiData.ubiScores(:, i);
    
    % Set COBRA solver parameters
    changeCobraSolverParams('LP', 'feasTol', 1e-9);

    % Extract the models
    extracted_models = extract_mCADRE_models(model, currentUbiScores, cell_line, phase);

    % Retrieve the extracted model
    r1 = extracted_models.red_models.rxns(any(extracted_models.red_models.retained_rxns, 2));
    reduced_model = removeRxns(model, model.rxns(~ismember(model.rxns, r1)));
    reduced_model = removeUnusedGenes(reduced_model);

    % Save the reduced model for this sample condition
    save(['reduced_model_' sample '.mat'], 'reduced_model');
end