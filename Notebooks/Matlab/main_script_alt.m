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


% Step 4: Define Protected Reactions List
protected_reactions = {
    'biomass_cho_s', 'DNAsyn', 'LipidSyn_cho_s', 'PROTsyn_cho_s', 'RNAsyn_cho_s', 'EX_bhb_e', 'EX_nh4_e', 'EX_ac_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asn_L_e', 'EX_asp_L_e', 'EX_2hb_e', 'EX_cit_e',
    'EX_cys_L_e', 'EX_etoh_e', 'EX_for_e', 'EX_fum_e', 'EX_glc_e', 'EX_glu_L_e', 'EX_gln_L_e', 'EX_gly_e', 'EX_his_L_e', 'EX_4hpro_e',
    'EX_ile_L_e', 'EX_lac_L_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_mal_L_e', 'EX_met_L_e', 'EX_phe_L_e', 'EX_pro_L_e', 'EX_5oxpro_e',
    'EX_pyr_e', 'EX_ser_L_e', 'EX_thr_L_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_val_L_e'
};

% Step 4: Model extraction using mCADRE
sampleConditions = UbiData.Condition;

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
    extracted_models = extract_mCADRE_models(model, currentUbiScores, cell_line, phase, protected_reactions);

    % Retrieve the extracted model
    r1 = extracted_models.red_models.rxns(any(extracted_models.red_models.retained_rxns, 2));
    reduced_model = removeRxns(model, model.rxns(~ismember(model.rxns, r1)));
    reduced_model = removeUnusedGenes(reduced_model);

    % Save the reduced model for this sample condition
    save(['reduced_model_' condition '.mat'], 'reduced_model');
end
