% Step 1: Initialize the COBRA toolbox
initCobraToolbox(0);
% 
% path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';
% addpath('/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells/Notebooks/Matlab/Model_Extraction');


path = '';
addpath('');

% Step 2: Ensure that the parent genome-scale model is loaded into the
% workspace and assigned the variable "model"
model = readCbModel(fullfile(path, 'Notebooks', 'iCHO3K.mat'));
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

confidenceScores = zeros(size(model.rxns));
UbiData.confidenceScores = confidenceScores;
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
protected_reactions = {
    'biomass_cho_s', 'DNAsyn', 'LipidSyn_cho_s', 'PROTsyn_cho_s', 'RNAsyn_cho_s', 'EX_bhb_e', 'EX_nh4_e', 'EX_ac_e', 'EX_ala_L_e', 'EX_arg_L_e', 'EX_asn_L_e', 'EX_asp_L_e', 'EX_2hb_e', 'EX_cit_e', ...
    'EX_cys_L_e', 'EX_etoh_e', 'EX_for_e', 'EX_fum_e', 'EX_glc_e', 'EX_glu_L_e', 'EX_gln_L_e', 'EX_gly_e', 'EX_his_L_e', 'EX_4hpro_e', ...
    'EX_ile_L_e', 'EX_lac_L_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_mal_L_e', 'EX_met_L_e', 'EX_phe_L_e', 'EX_pro_L_e', 'EX_5oxpro_e', 'EX_pyr_e', ...
    'EX_ser_L_e', 'EX_thr_L_e', 'EX_trp_L_e', 'EX_tyr_L_e', 'EX_val_L_e', 'EX_h2o_e', 'EX_h_e', 'EX_o2_e', 'EX_hco3_e', 'EX_so4_e', 'EX_pi_e', ...
    'SK_Asn_X_Ser_Thr_r', 'SK_Tyr_ggn_c', 'SK_Ser_Thr_g', 'SK_pre_prot_r'
};


% Step 5: Load the uptake and secretion files
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
        extracted_models = extract_mCADRE_models(model, currentUbiScores, cell_line, phase, protected_reactions, 0);
%         extracted_models = extract_mCADRE_models(model, currentUbiScores, cell_line, phase, protected_reactions, 1);
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

%% Example code
% intervals = {'P0', 'P2'; 'P2', 'P4'; 'P4', 'P6'; 'P6', 'P8'; 'P8', 'P12'; 'P12', 'P14'};
% preprocessing = 'raw'; % 'smoothened'
% 
% % Ensure results structures are initialized
% results = [];
% results_for_check = [];
% cumulative_results = containers.Map;
% 
% 
% % List of batch IDs
% batch_ids = {'U1', 'U2', 'U3', 'U4', 'U5', 'U6', 'U7', 'U8'};
% % Initialize result table
% results = table();
% results_for_check = table();
% 
% for i = 1:length(batch_ids)
%     batch_id = batch_ids{i};
%     if isKey(batch_dfs_meta, batch_id)
%         % Extract the table from the map
%         batch_data = batch_dfs_meta(batch_id);
%         
%         % Drop the first row
%         batch_data(1, :) = [];
%         
%         % Update the map with the modified table
%         batch_dfs_meta(batch_id) = batch_data;
%         
%         % Get the list of metabolites
%         met_list = batch_data.Properties.VariableNames(3:end);
% 
%         % Calculate the specific rates for each interval
%         for j = 1:size(intervals, 1)
%             start = intervals{j, 1};
%             end_ = intervals{j, 2};
%             interval_key = [start '_' end_];  % Create a unique key for the interval
% 
%             end_row = batch_data(strcmp(batch_data.SampleID, end_), :);
%             start_row = batch_data(strcmp(batch_data.SampleID, start), :);
% 
%             % Proceed only if both rows exist
%             if ~isempty(end_row) && ~isempty(start_row)
%                 try
%                     if strcmp(preprocessing, 'smoothened')
%                         % Extract necessary values from the smoothened data
%                         batch_data_smooth = batch_dfs(batch_id); % Assuming batch_dfs is a map with smoothened data
%                         growth_rate = batch_data_smooth(strcmp(batch_data_smooth.SampleID, end_), :).GrowthRate;
%                         VC_t = batch_data_smooth(strcmp(batch_data_smooth.SampleID, end_), :).fitted_VCD;
%                         VC_i = batch_data_smooth(strcmp(batch_data_smooth.SampleID, start), :).fitted_VCD;
%                         IVCD = [];
% 
%                         current_time = batch_data_smooth(strcmp(batch_data_smooth.SampleID, end_), :).Age_h_;
%                         initial_time = batch_data_smooth(strcmp(batch_data_smooth.SampleID, start), :).Age_h_;
%                     elseif strcmp(preprocessing, 'raw')
%                         % Extract necessary values from the raw data
%                         batch_data_raw = batch_dfs_raw(batch_id); % Assuming batch_dfs_raw is a map with raw data
%                         growth_rate = batch_data_raw(strcmp(batch_data_raw.SampleID, end_), :).GrowthRate;
%                         IVCD = batch_data_raw(strcmp(batch_data_raw.SampleID, end_), :).VCD;
%                         VC_t = [];
%                         VC_i = [];
% 
%                         current_time = batch_data_raw(strcmp(batch_data_raw.SampleID, end_), :).Age_h_;
%                         initial_time = batch_data_raw(strcmp(batch_data_raw.SampleID, start), :).Age_h_;
%                     end
% 
%                     % Save the growth rate in the cumulative results
%                     new_entry = table({batch_id}, {interval_key}, {'SGR'}, growth_rate, ...
%                                       'VariableNames', {'BatchID', 'Interval', 'Metabolite', 'SpecificRate'});
%                     results = [results; new_entry];
% 
%                     % Loop through each metabolite to calculate specific rates
%                     for k = 1:length(met_list)
%                         met_name = met_list{k};
%                         delta_met = end_row.(met_name) - start_row.(met_name);
% 
%                         % Calculate specific rate
%                         specific_rate = calculate_specific_rate(delta_met, growth_rate, initial_time, current_time, VC_t, VC_i, IVCD);
% 
%                         % Save the results
%                         new_entry = table({batch_id}, {interval_key}, {met_name}, specific_rate, ...
%                                           'VariableNames', {'BatchID', 'Interval', 'Metabolite', 'SpecificRate'});
%                         results = [results; new_entry];
% 
%                         % Save the check results
%                         check_entry = table({num2str(delta_met)}, {num2str(growth_rate)}, {num2str(VC_t)}, ...
%                                             {num2str(VC_i)}, {num2str(IVCD)}, {num2str(initial_time)}, {num2str(current_time)}, ...
%                                             'VariableNames', {'DeltaMet', 'GrowthRate', 'VC_t', 'VC_i', 'IVCD', 'InitialTime', 'CurrentTime'});
%                         results_for_check = [results_for_check; check_entry];
%                     end
%                 catch ME
%                     fprintf('Error processing interval %s-%s for batch %s: %s\n', start, end_, batch_id, ME.message);
%                 end
%             end
%         end
%     end
% end
% 
% % Display the cumulative results for verification
% disp(results);
% disp(results_for_check);
