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
    'S1', 'WT_P2_Bio141'; 'S7', 'WT_P2_Bio142'; 'S13', 'WT_P2_Bio143';
    'S2', 'WT_P4_Bio141'; 'S8', 'WT_P4_Bio142'; 'S14', 'WT_P4_Bio143';
    'S3', 'WT_P6_Bio141'; 'S9', 'WT_P6_Bio142'; 'S15', 'WT_P6_Bio143';
    'S4', 'WT_P8_Bio141'; 'S10', 'WT_P8_Bio142'; 'S16', 'WT_P8_Bio143';
    'S5', 'WT_P12_Bio141';'S11', 'WT_P12_Bio142';
    'S6', 'WT_P14_Bio141';'S12', 'WT_P14_Bio142';'S17', 'WT_P14_Bio143';
    'S18', 'ZeLa_P4_Bio144'; 'S23', 'ZeLa_P4_Bio145'; 'S28', 'ZeLa_P4_Bio146'; 'S34', 'ZeLa_P4_Bio147'; 'S39', 'ZeLa_P4_Bio148';
    'S19', 'ZeLa_P6_Bio144'; 'S29', 'ZeLa_P6_Bio146'; 'S35', 'ZeLa_P6_Bio17'; 'S40', 'ZeLa_P6_Bio148';
    'S20', 'ZeLa_P8_Bio144'; 'S24', 'ZeLa_P8_Bio145'; 'S30', 'ZeLa_P8_Bio146'; 'S36', 'ZeLa_P8_Bio147'; 'S41', 'ZeLa_P8_Bio148';
    'S25', 'ZeLa_P12_Bio145';'S31', 'ZeLa_P12_Bio146';'S42', 'ZeLa_P12_Bio148';
    'S21', 'ZeLa_P14_Bio144';'S26', 'ZeLa_P14_Bio145';'S32', 'ZeLa_P14_Bio146';'S37', 'ZeLa_P14_Bio147';'S43', 'ZeLa_P14_Bio148';
    'S22', 'ZeLa_P2_Bio145'; 'S27', 'ZeLa_P2_Bio146'; 'S33', 'ZeLa_P2_Bio147'; 'S38', 'ZeLa_P2_Bio148';
};


for i = 1:length(UbiData.ubiScores(1,:))
    % Determine the corresponding cell line and phase
    sample = sampleConditions{i, 1};
    condition = sampleConditions{i, 2};
    splitCondition = strsplit(condition, '_');
    cell_line = splitCondition{1};
    phase = [splitCondition{2} '_' splitCondition{3}];

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
