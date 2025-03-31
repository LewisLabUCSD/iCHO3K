% Step 1: Initialize the COBRA toolbox
initCobraToolbox(0);

% path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';
% addpath('/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells/Notebooks/Matlab/Model_Extraction');
path = 'C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells';
addpath('C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells\Notebooks\Matlab\Model_Extraction');

% Step 2: Ensure that the parent genome-scale model is loaded into the
% workspace and assigned the variable "model"
model = readCbModel(fullfile(path, 'Notebooks', 'iCHO3K_unblocked.mat'));

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
% Option1: Use zero confidence scores
confidenceScores = zeros(size(model.rxns));
UbiData.confidenceScores = confidenceScores;

% Option2: Use standep-generated confidence scores
% % Check if the confidence score exists in the loaded data
% if isfield(UbiData, 'confidenceScores')
%     confidenceScores = UbiData.confidenceScores;    
% else
%     % Load confidence scores from CSV file
%     csvFilePath = fullfile(path, 'Data/Context_specific_models', 'confidence_scores.csv');
%     csvData = readtable(csvFilePath);
%     
%     % Ensure rxn exists in the loadedData
%     if isfield(UbiData, 'rxns')
%         rxns = UbiData.rxns;
%         
%         % Check if the lengths of rxn and confidenceScores are equal
%         if height(csvData) == length(rxns)
%             UbiData.confidenceScores = csvData.Var1; % Assuming the scores are in the first column
%             fprintf('The confidence scores were successfully added to loadedData.\n');
%         else
%             error('The length of "confidenceScores" from CSV does not match the length of "rxn".');
%         end
%     else
%         error('The variable "rxns" does not exist in the UbiData.');
%     end
% end

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
% protected_reactions = {
%     'biomass_cho_s', 'DNAsyn', 'LipidSyn_cho_s', 'PROTsyn_cho_s', 'RNAsyn_cho_s', 'EX_bhb_e', 'EX_nh4_e', 'EX_ac_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asn_L_e', 'EX_asp_L_e', 'EX_2hb_e', 'EX_cit_e', ...
%     'EX_cys_L_e', 'EX_etoh_e', 'EX_for_e', 'EX_fum_e', 'EX_glc_e', 'EX_glu_L_e', 'EX_gln_L_e', 'EX_gly_e', 'EX_his_L_e', 'EX_4hpro_e', ...
%     'EX_ile_L_e', 'EX_lac_L_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_mal_L_e', 'EX_met_L_e', 'EX_phe_L_e', 'EX_pro_L_e', 'EX_5oxpro_e', 'EX_pyr_e', ...
%     'EX_ser_L_e', 'EX_thr_L_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_val_L_e', 'EX_h2o_e', 'EX_h_e', 'EX_o2_e', 'EX_hco3_e', 'EX_so4_e', 'EX_pi_e', ...
%     'SK_Asn_X_Ser_Thr_r', 'SK_Tyr_ggn_c', 'SK_Ser_Thr_g', 'SK_pre_prot_r'
% };
% 
% protected_reactions = {
%     'biomass_cho_s', 'DNAsyn', 'LipidSyn_cho_s', 'PROTsyn_cho_s', 'RNAsyn_cho_s', 'EX_bhb_e', 'EX_nh4_e', 'EX_ac_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asn_L_e', 'EX_asp_L_e', 'EX_2hb_e', 'EX_cit_e', ...
%     'EX_cys_L_e', 'EX_etoh_e', 'EX_for_e', 'EX_fum_e', 'EX_glc_e', 'EX_glu_L_e', 'EX_gln_L_e', 'EX_gly_e', 'EX_his_L_e', 'EX_4hpro_e', ...
%     'EX_ile_L_e', 'EX_lac_L_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_mal_L_e', 'EX_met_L_e', 'EX_phe_L_e', 'EX_pro_L_e', 'EX_5oxpro_e', 'EX_pyr_e', ...
%     'EX_ser_L_e', 'EX_thr_L_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_val_L_e', 'EX_h2o_e', 'EX_h_e', 'EX_o2_e', 'EX_hco3_e', 'EX_so4_e', 'EX_pi_e'
% };

protected_reactions = {
'biomass_cho_s', 'DNAsyn', 'LipidSyn_cho_s', 'PROTsyn_cho_s', 'RNAsyn_cho_s', ...
'EX_bhb_e', 'EX_nh4_e', 'EX_ac_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asn_L_e', ...
'EX_asp_L_e', 'EX_2hb_e', 'EX_cit_e', 'EX_cys_L_e', 'EX_etoh_e', 'EX_for_e', ...
'EX_fum_e', 'EX_glc_e', 'EX_glu_L_e', 'EX_gln_L_e', 'EX_gly_e', 'EX_his_L_e', ...
'EX_4hpro_e', 'EX_ile_L_e', 'EX_lac_L_e', 'EX_leu_L_e', 'EX_lys_L_e', ...
'EX_mal_L_e', 'EX_met_L_e', 'EX_phe_L_e', 'EX_pro_L_e', 'EX_5oxpro_e', ...
'EX_pyr_e', 'EX_ser_L_e', 'EX_thr_L_e', 'EX_trp_L_e', 'EX_tyr_L_e', ...
'EX_val_L_e', 'EX_h2o_e', 'EX_h_e', 'EX_o2_e', 'EX_hco3_e', 'EX_so4_e', ...
'EX_pi_e', 'SK_Asn_X_Ser_Thr_r', 'SK_Tyr_ggn_c', 'SK_Ser_Thr_g', 'SK_pre_prot_r', ...
'ACONTa', 'ACONTam', 'ACONTb', 'ACONTbm', 'ACITL', 'ACONT', 'ACONTm', 'AKGDm', ...
'CITL', 'CSm', 'FUM', 'FUMm', 'ICDHxm', 'ICDHy', 'ICDHyp', 'ICDHyrm', 'MDH', ...
'MDHm', 'MDHx', 'r0081', 'SUCD1m', 'SUCOAS1m', 'SUCOASm'
};

% Step 5: Load the Python dictionary for reaction bounds
% Python code to convert the pickle file to a .mat file:
% ------------------------------------------------------------------------------------
% import pickle
% import scipy.io
% 
% # Load the pickle file
% with open(r'path and name for your pkl file', 'rb') as f:
%     data = pickle.load(f)
% 
% # Save the data to a .mat file for MATLAB
% scipy.io.savemat('uptake_secretion_wt.mat', {'data': data})
% ------------------------------------------------------------------------------------
%
% Load the experimental data dictionaries
uptsec_intrvl_wt = load('uptake_secretion_wt.mat') ;
uptsec_intrvl_zela = load('uptake_secretion_zela.mat');

% Step 6: Constraining the reference models with experimental (simulated)
% datasets
% Initialize a structure to hold all models
all_models = struct();

% Loop through ubiquity scores
for i = 1:length(UbiData.ubiScores(1,:))
    sampleConditions = UbiData.Condition;
    condition = sampleConditions{1, i};
    fprintf('Successfully loaded condition "%s".\n', condition);
    
    splitCondition = strsplit(condition, '_');
    cell_line = splitCondition{1};
    phase = [splitCondition{2} '_' splitCondition{3}];
    
    % Select the appropriate dictionary based on the cell line name
    if startsWith(cell_line, 'WT')
        bounds_dict = uptsec_intrvl_wt;
    elseif strcmp(cell_line, 'ZeLa')
        bounds_dict = uptsec_intrvl_zela;
    else
        error('Unknown cell line: %s', cell_line);
    end
    
    % Determine the time based on the phase
    phase_start = phase(1:3);
    switch phase_start
        case 'P2_'
            time = 'P0 to P2';
        case 'P4_'
            time = 'P2 to P4';
        case 'P6_'
            time = 'P4 to P6';
        case 'P8_'
            time = 'P6 to P8';
        case 'P12'
            time = 'P8 to P12';
        case 'P14'
            time = 'P12 to P14';
        otherwise
            error('Unknown phase: %s', phase);
    end
    
    % Update the ubiquity scores for the current sample
    currentUbiScores = UbiData;
    currentUbiScores.ubiScores = UbiData.ubiScores(:, i);

    % Constrain the bounds of reactions based on the experimental data
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
            % Add any additional handling if necessary
        end
        sink_reactions = {'SK_Asn_X_Ser_Thr_r', 'SK_Tyr_ggn_c', 'SK_Ser_Thr_g', 'SK_pre_prot_r'};
        for sk = 1:length(sink_reactions)
            if strcmp(rxn, sink_reactions{sk}) % Set constrains for sink reactions
                model = changeRxnBounds(model, rxn, -0.1, 'l');  % Lower bound
                model = changeRxnBounds(model, rxn, 1000, 'u');  % Upper bound
            end
        end
        
    end

    % Add cell_line and phase fields to the model
    model.cell_line = cell_line;
    model.phase = phase;
    
    % Save the model to a structure with a unique name
    model_name = sprintf('%s_%s_model', cell_line, phase);
    all_models.(model_name) = model;
end
% Save the structure to a .mat file
save('all_models_NOV.mat', 'all_models');

% Step 7: Load the constrained models and identify non-contextualing
% models with NaN flux of objective function
% Load the saved structure to verify
loaded_data = load('all_models_NOV.mat');
loaded_models = loaded_data.all_models;

% Display the names of the models saved
model_names = fieldnames(loaded_models);
fprintf('The following models were saved:\n');

% Initialize a cell array to hold models with NaN objective values
nan_models = {};
for i = 1:length(model_names)
    fprintf('%s\n', model_names{i});
end

% Verify the bounds of each reaction in one of the models
for i = 1:length(model_names)
    example_model_name = model_names{i};
    example_model = loaded_models.(example_model_name);
    
    % Display some details of the example model
    fprintf('Details of model %s:\n', example_model_name);
    
    % Verify reaction bounds for biomass, etoh, and bounds_dict reactions
    reactions_to_check = {'biomass_cho_s', 'EX_etoh_e'};
    for j = 1:length(reactions_to_check)
        rxn = reactions_to_check{j};
        if any(strcmp(example_model.rxns, rxn))
            idx = find(strcmp(example_model.rxns, rxn));
            lb = example_model.lb(idx);
            ub = example_model.ub(idx);
            fprintf('Reaction %s: lower bound = %f, upper bound = %f\n', rxn, lb, ub);
        end
    end
    
    % Check reactions in bounds_dict
    bounds_dict_reactions = fieldnames(bounds_dict);
    for j = 1:length(bounds_dict_reactions)
        rxn = bounds_dict_reactions{j};
        if any(strcmp(example_model.rxns, rxn))
            idx = find(strcmp(example_model.rxns, rxn));
            lb = example_model.lb(idx);
            ub = example_model.ub(idx);
            fprintf('Reaction %s: lower bound = %f, upper bound = %f\n', rxn, lb, ub);
        end
    end
    % Set biomass as the objective and optimize the model
    example_model = changeObjective(example_model, 'biomass_cho_s');
    sol = optimizeCbModel(example_model);
    fprintf('Optimization result for model %s: Objective value = %f\n', example_model_name, sol.f);
    
        % Check if the objective value is NaN and save the model if it is
    if isnan(sol.f)
        nan_models{i,1} = example_model; % Save the model in the cell array
        fprintf('Model %s has NaN objective value and was added to nan_models.\n', example_model_name);
    end
    
end


% Save nan_models to a .mat file for later use
save('nan_models_NOV.mat', 'nan_models');

% Display the number of models with NaN objective values
fprintf('Total number of models with NaN objective values: %d\n', length(nan_models));

% Step 8: Model extraction using mCADRE -- We can start with all the P4 models (P2 to P4). And then we can move to P2 models (P0 to P2)
% for i = 1:length(UbiData.ubiScores(1,:))
for i = 1:length(UbiData.ubiScores(1,:))
    % Determine the corresponding cell line and phase
    condition = sampleConditions{1, i};
    fprintf('Successfully loaded condition "%s".\n', condition);
    splitCondition = strsplit(condition, '_');
    cell_line = splitCondition{1};
    phase = [splitCondition{2}, '_', splitCondition{3}];

    % Add filtering on P4 cell lines
    if startsWith(phase, 'P4') || startsWith(phase, 'P6')
        % Load the appropriate model from all_models
        model_name = sprintf('%s_%s_model', cell_line, phase);
        CBmodel = loaded_models.(model_name);
        
        % Set COBRA solver parameters
        changeCobraSolverParams('LP', 'feasTol', 1e-9);
        currentUbiScores.ubiScores = UbiData.ubiScores(:, i);   
        % Extract the models
        extracted_models = extract_mCADRE_models_v1_1(CBmodel, currentUbiScores, cell_line, phase, protected_reactions, 0);
        % Retrieve the extracted model
        r1 = extracted_models.red_models.rxns(any(extracted_models.red_models.retained_rxns, 2));
        reduced_model = removeRxns(CBmodel, CBmodel.rxns(~ismember(CBmodel.rxns, r1)));
        reduced_model = removeUnusedGenes(reduced_model);

        % Save the reduced model for this sample condition
        save(['NOV_reduced_model_CF_' condition '.mat'], 'reduced_model');
    elseif (startsWith(phase, 'P2') || startsWith(phase, 'P8')) && startsWith(cell_line, 'ZeLa')
        % Load the appropriate model from all_models
        model_name = sprintf('%s_%s_model', cell_line, phase);
        CBmodel = loaded_models.(model_name);
        
        % Set COBRA solver parameters
        changeCobraSolverParams('LP', 'feasTol', 1e-9);

        % Extract the models
        extracted_models = extract_mCADRE_models_v1_1(CBmodel, currentUbiScores, cell_line, phase, protected_reactions, 0);
        % Retrieve the extracted model
        r1 = extracted_models.red_models.rxns(any(extracted_models.red_models.retained_rxns, 2));
        reduced_model = removeRxns(CBmodel, CBmodel.rxns(~ismember(CBmodel.rxns, r1)));
        reduced_model = removeUnusedGenes(reduced_model);
        reduced_model.zeroExpRxns = extracted_models.red_models.zeroExpRxns;

        % Save the reduced model for this sample condition
        save(['NOV_reduced_model_CF_' condition '.mat'], 'reduced_model');
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

% 
% % initCobraToolbox(0);
% load('AUG15_reduced_model_CF_WT_P4_Bio143.mat')
% % Step 1: Optimize the model
% ans = optimizeCbModel(reduced_model, 'max', 'one');
% % Step 2: Extract relevant information
% reactionIDs = reduced_model.rxns; % Reaction IDs
% reactionFormulas = printRxnFormula(reduced_model); % Reaction formulas
% reactionSubsystems = reduced_model.subSystems; % Reaction subsystems
% fluxValues = ans.x; % Flux values from the optimization result
% % Step 3: Combine the extracted information into a table
% resultTable = table(reactionIDs, reactionFormulas, reactionSubsystems, fluxValues, ...
% 'VariableNames', {'ReactionID', 'ReactionFormula', 'Subsystem', 'FluxValue'});
% % Step 4: Export the table to an Excel file
% writetable(resultTable, 'optimization_results_0820.xlsx');
% 
% 
% for i = 1:length(protected_reactions)
% reaction_to_check = protected_reactions{i};
% if ismember(reaction_to_check, reduced_model.rxns)
%     fprintf('The reaction "%s" is present in reduced_model.rxns.\n', reaction_to_check);
% else
%     fprintf('The reaction "%s" is not present in reduced_model.rxns.\n', reaction_to_check);
% end
% end