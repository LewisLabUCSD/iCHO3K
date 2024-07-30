% Step 1: Initialize the COBRA toolbox
% initCobraToolbox(0);
% 
% path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';
% addpath('/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells/Notebooks/Matlab/Model_Extraction');


path = 'C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells';
addpath('C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells\Notebooks\Matlab\Model_Extraction');

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


% Check if the confidence score exists in the loaded data
if isfield(UbiData, 'confidenceScores')
    confidenceScores = UbiData.confidenceScores;    
else
    % Load confidence scores from CSV file
    csvFilePath = fullfile(path, 'Data/Context_specific_models', 'confidence_scores.csv');
    csvData = readtable(csvFilePath);
    
    % Ensure rxn exists in the loadedData
    if isfield(UbiData, 'rxns')
        rxns = UbiData.rxns;
        
        % Check if the lengths of rxn and confidenceScores are equal
        if height(csvData) == length(rxns)
            UbiData.confidenceScores = csvData.Var1; % Assuming the scores are in the first column
            fprintf('The confidence scores were successfully added to loadedData.\n');
        else
            error('The length of "confidenceScores" from CSV does not match the length of "rxn".');
        end
    else
        error('The variable "rxns" does not exist in the UbiData.');
    end
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
    'biomass_cho_s', 'DNAsyn', 'LipidSyn_cho_s', 'PROTsyn_cho_s', 'RNAsyn_cho_s', 'EX_bhb_e', 'EX_nh4_e', 'EX_ac_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asn_L_e', 'EX_asp_L_e', 'EX_2hb_e', 'EX_cit_e', ...
    'EX_cys_L_e', 'EX_etoh_e', 'EX_for_e', 'EX_fum_e', 'EX_glc_e', 'EX_glu_L_e', 'EX_gln_L_e', 'EX_gly_e', 'EX_his_L_e', 'EX_4hpro_e', ...
    'EX_ile_L_e', 'EX_lac_L_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_mal_L_e', 'EX_met_L_e', 'EX_phe_L_e', 'EX_pro_L_e', 'EX_5oxpro_e', 'EX_pyr_e', ...
    'EX_ser_L_e', 'EX_thr_L_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_val_L_e', 'EX_h2o_e', 'EX_h_e', 'EX_o2_e', 'EX_hco3_e', 'EX_so4_e', 'EX_pi_e', ...
    'SK_Asn_X_Ser_Thr_r', 'SK_Tyr_ggn_c', 'SK_Ser_Thr_g', 'SK_pre_prot_r'
};


% Step 5: Load the Python dictionary for reaction bounds
% % Set up Python environment (using Anaconda Python environment)
% pyenv('Version', 'C:\Users\user\anaconda3\python.exe');
% 
% % Import Python modules
% pickle = py.importlib.import_module('pickle');
% builtins = py.importlib.import_module('builtins');
% 
% % Load pickle files
% uptsec_intrvl_wt = pickle.load(builtins.open('..\Data\Uptake_Secretion_Rates\uptake_secretion_intrvl_wt_dict.pkl', 'rb'));
% uptsec_intrvl_zela = pickle.load(builtins.open('..\Data\Uptake_Secretion_Rates\uptake_secretion_intrvl_zela_dict.pkl', 'rb'));
% 
% % Display the loaded data
% disp(uptsec_intrvl_wt);
% disp(uptsec_intrvl_zela);

uptsec_intrvl_wt = load('uptake_secretion_wt.mat');
uptsec_intrvl_zela = load('uptake_secretion_zela.mat');


% Step 6: Model extraction using mCADRE -- We can start with all the P4 models (P2 to P4). And then we can move to P2 models (P0 to P2)
sampleConditions = UbiData.Condition;

% for i = 1:length(UbiData.ubiScores(1,:))
for i = 1:length(UbiData.ubiScores(1,:))
    % Determine the corresponding cell line and phase
    condition = sampleConditions{1, i};
    fprintf('Successfully loaded condition "%s".\n', condition);
    splitCondition = strsplit(condition, '_');
    cell_line = splitCondition{1};
    phase = [splitCondition{2} '_' splitCondition{3}];

    % Add filtering on P4 cell lines
    if startsWith(phase, 'P4')
        % Select the appropriate dictionary based on the cell line name
        if startsWith(cell_line, 'WT')
            bounds_dict = uptsec_intrvl_wt;
        elseif strcmp(cell_line, 'ZeLa')
            bounds_dict = uptsec_intrvl_zela;
        else
            error('Unknown cell line: %s', cell_line);
        end
        fprintf('Successfully loaded phase for "%s".\n', phase);
        % Update the ubiquity scores for the current sample
        currentUbiScores = UbiData;
        currentUbiScores.ubiScores = UbiData.ubiScores(:, i);

        % Constrain the bounds of reactions based on the experimental data
        time = 'P2 to P4';
        for j = 1:length(model.rxns)
            rxn = model.rxns{j};
            if strcmp(rxn, 'biomass_cho_s') % Set constrains for biomass
                bounds = bounds_dict.exp_growth_rate.(time);
                model = changeRxnBounds(model, rxn, bounds(1), 'l');  % Lower bound
                model = changeRxnBounds(model, rxn, bounds(2), 'u');  % Upper bound
            end
            if strcmp(rxn, 'EX_etoh_e') % Models are not feasible when forced to secrete ethanol
                model = changeRxnBounds(model, rxn, -0.1, 'l');  % Lower bound
                model = changeRxnBounds(model, rxn, 0.1, 'u');  % Upper bound
           elseif isfield(bounds_dict, rxn)
                bounds = bounds_dict.(rxn).(time);
                model = changeRxnBounds(model, rxn, bounds(1), 'l');  % Lower bound
                model = changeRxnBounds(model, rxn, bounds(2), 'u');  % Upper bound
            else
                
            end
        end

        % Set COBRA solver parameters
        changeCobraSolverParams('LP', 'feasTol', 1e-9);


        % Extract the models
%         extracted_models = extract_mCADRE_models(model, currentUbiScores, cell_line, phase, protected_reactions, 0);
        extracted_models = extract_mCADRE_models(model, currentUbiScores, cell_line, phase, protected_reactions, 1);
        % Retrieve the extracted model
        r1 = extracted_models.red_models.rxns(any(extracted_models.red_models.retained_rxns, 2));
        reduced_model = removeRxns(model, model.rxns(~ismember(model.rxns, r1)));
        reduced_model = removeUnusedGenes(reduced_model);

        % Save the reduced model for this sample condition
        save(['reduced_model_CF_' condition '.mat'], 'reduced_model');
    end
end


for i = 1:length(protected_reactions)
reaction_to_check = protected_reactions{i};
if ismember(reaction_to_check, reduced_model.rxns)
    fprintf('The reaction "%s" is present in reduced_model.rxns.\n', reaction_to_check);
else
    fprintf('The reaction "%s" is not present in reduced_model.rxns.\n', reaction_to_check);
end
end