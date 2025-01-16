% Load the experimental data dictionaries
uptsec_intrvl_wt = load('uptake_secretion_wt.mat');
uptsec_intrvl_zela = load('uptake_secretion_zela.mat');

% Initialize a structure to hold all models
all_models = struct();

% Specify the directory containing your context-specific models
modelDir = 'C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells\Data\Context_specific_models\ecModels';  % Update this to the correct path
keyword = 'P6';
% Get all .mat files in the directory
modelFiles = dir(fullfile(modelDir, '*.mat'));
% Filter files to include only those containing the keyword in their name
modelFiles = modelFiles(contains({modelFiles.name}, keyword));
ref_model = load('iCHO3K.mat');
identifier = "New";
% removeReactions = ["OCOAT1m"];
removeReactions_ATP = ["ICDHyr", "OCOAT1m", "r0082", "r0083", "SCP22x", "TMNDNCCOAtx", "OCCOAtx", "r0391", "BiGGRxn67", "r2247", "r2246", "r2245", "r2317", "HMR_0293", "HMR_7741", "r1453", "r1393", "NICRNS", "GapFill-R08726", "RE2915M", "HMR_3288", "HMR_1325", "RE2439C", "r1450", "RE3477C", "AAPSAS", "r0698", "3HDH260p", "HMR_3272", "ACOAD183n3m", "GapFill-R01463", "r1468", "r0655", "r0603", "r0541", "HMR_1329", "GapFill-R03599", "OIVD1m", "OIVD3m", "2OXOADOXm", "r0386", "r0451", "GLYCLm", "MACACI", "r0509", "r0425", "r0423", "r0424"];

addReactionId = ["AKGDm"];
checkReactions = ["ACONTa", "ACONTam", "ACONTb", "ACONTbm", "ACITL", ...
    "ACONT", "ACONTm", "AKGDm", "CITL", "CSm", "FUM", "FUMm", ...
    "ICDHxm", "ICDHy", "ICDHyp", "ICDHyrm", "MDH", "MDHm", "MDHx",... 
    "r0081", "SUCD1m", "SUCOAS1m", "SUCOASm"...
];

% Loop through each model file
for i = 1:length(modelFiles)
    % Load the model file
    modelFilePath = fullfile(modelDir, modelFiles(i).name);
    loaded_model = load(modelFilePath);
    % Assuming the loaded model is stored in a variable called 'model'
    if isfield(loaded_model, 'model')
        model = loaded_model.model;
        
        % Extract metadata from the filename if necessary
        [~, fileName, ~] = fileparts(modelFiles(i).name);  % Get the filename without extension
        fprintf('Processing model: %s\n', fileName);
        % Define sampling options
        options.numSamples = 5000;  % Adjust the number of samples as needed
        options.stepsPerPoint = 100;
        options.populationScale = 2;
%         model_1 = model
        model_1 = changeModelCondition(model, [1,2,3], removeReactions_ATP, checkReactions, checkReactions, ref_model.iCHO3K);
        model_1 = changeObjective(model_1, 'biomass_cho_s');
%         model_1 = changeObjective(model_1, 'AKGDm');
        
        %% Optional!!
        model_1.lb(strcmp(model.rxns, 'ATPM')) = 0; 
        model_1.ub(strcmp(model.rxns, 'ATPM')) = 1000;
%%
        % Run the sampling for the loaded model
%         sample = fluxSamplerADSB(model_1, options);
        sample = fluxSamplerADSB_v2(model_1, options);
        fprintf('\nChecking specified reactions:\n');
        for i3 = 1:length(checkReactions)
            reaction = checkReactions{i3};
            if ismember(reaction, sample.rxns)
                fprintf('Reaction %s is present in the model.\n', reaction);
            else
                fprintf('Reaction %s is NOT present in the model.\n', reaction);
            end
        end
        for i4 = 1:length(removeReactions_ATP)
            reaction = removeReactions_ATP{i4};
            if ismember(reaction, sample.rxns)
                fprintf('Reaction %s is present in the model.\n', reaction);
            else
                fprintf('Reaction %s is NOT present in the model.\n', reaction);
            end
        end
        % Save the sampling result to a new .mat file with model context
        save_filename = fullfile('C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells\Data\Context_specific_models\sampling_ADSB_tmp', sprintf('ADSB_v2_sample_%s_%s.mat',identifier, fileName));
        save(save_filename, 'sample', 'fileName');
        fprintf('Saved sampling result to %s\n', save_filename);
        
        % After sampling is completed
        mean_sample = mean(sample.points, 2);  % Calculate mean of each reaction's flux across samples
        std_sample = std(sample.points, 0, 2);  % Calculate standard deviation for each reaction's flux

        % Get reaction names and formulas from model_1
        reaction_names = model_1.rxnNames;  % Reaction names
        reaction_formulas = printRxnFormula(model_1, 'rxnAbbrList', sample.rxns, 'printFlag', false);  % Reaction formulas

        % Create a table to store reaction information, including names, formulas, mean, and std
        reaction_summary = table(sample.rxns, reaction_names, reaction_formulas, mean_sample, std_sample, ...
                                 'VariableNames', {'ReactionID', 'ReactionName', 'Formula', 'Mean_Flux', 'Std_Flux'});

        % Save the summary table to an Excel file
        output_filename = fullfile('C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\temp_code\sampled_results', sprintf('ADSB_v2_reaction_summary_%s_%s.xlsx', identifier, fileName));
        writetable(reaction_summary, output_filename);

        fprintf('Summary data exported to %s\n', output_filename);
    else
        fprintf('No model found in %s\n', modelFiles(i).name);
    end
end

function model_1 = changeModelCondition(model, option, removeReactions, addReactionId, checkReactions, ref_model)
    model_1 = model;
    % Option 1: Remove specific reactions
%     if ismember(1, option)
%         for i = 1:length(removeReactions)
%             reactionId = removeReactions{i};
%             reactionIndex = find(strcmp(model_1.rxns, reactionId),1);
%             if ~isempty(reactionIndex)
%                 model_1 = removeRxns(model_1, reactionId);
%                 fprintf('Removed reaction: %s\n', reactionId);
%             else
%                 fprintf('Reaction %s not found.\n', reactionId);
%             end
%         end
%     end
    
    if ismember(1, option)
        for i = 1:length(removeReactions)
            reactionId = removeReactions{i};
            reactionIndex = find(strcmp(model_1.rxns, reactionId), 1);
            if ~isempty(reactionIndex)
                % Create a temporary copy of the model (to preserve the original model)
                temp_model = model_1;

                % Attempt to remove the reaction
                temp_model = removeRxns(temp_model, reactionId);
                fprintf('Attempting to remove reaction: %s\n', reactionId);

                % Optimize the temporary model after removing the reaction
                solution = optimizeCbModel(temp_model);

                % Check optimization result
                if solution.stat == 1 && any(solution.x ~= 0)
                    % Apply removal only if the solution is feasible and flux is not zero
                    model_1 = temp_model;
                    fprintf('Removed reaction: %s\n', reactionId);
                else
                    % Skip removal if solution is infeasible or all fluxes are zero
                    if solution.stat ~= 1
                        fprintf('Removal of reaction %s resulted in an infeasible solution. Skipping removal.\n', reactionId);
                    else
                        fprintf('Removal of reaction %s resulted in a zero flux solution. Skipping removal.\n', reactionId);
                    end
                end
            else
                fprintf('Reaction %s not found.\n', reactionId);
            end
        end
    end


    % Option 2: Add specific reactions from ref_model
    if ismember(2, option)
        % Loop through each reaction in addReactionId
        for i2 = 1:length(addReactionId)
            rxn_add_id = addReactionId{i2};

            % Check if the reaction already exists in model_1
            if ~ismember(rxn_add_id, model_1.rxns)
                fprintf('Reaction %s not found in model_1. Attempting to add it from ref_model...\n', rxn_add_id);

                % Find the reaction in ref_model
                refIndex = find(strcmp(ref_model.rxns, rxn_add_id), 1);
                if ~isempty(refIndex)
                    fprintf('Reaction %s found in ref_model. Proceeding with addition...\n', rxn_add_id);

                    % Extract metabolites and stoichiometry for the reaction
                    metIndices = find(ref_model.S(:, refIndex) ~= 0);  % Indices of non-zero stoichiometric coefficients
                    metaboliteList = ref_model.mets(metIndices);  % Corresponding metabolites
                    stoichCoeffList = full(ref_model.S(metIndices, refIndex));  % Corresponding stoichiometric coefficients

                    % Use COBRA Toolbox's addReaction function to add the reaction to model_1
                    model_1 = addReaction(model_1, ...
                                          rxn_add_id, ...
                                          'metaboliteList', metaboliteList, ...
                                          'stoichCoeffList', stoichCoeffList, ...
                                          'lowerBound', ref_model.lb(refIndex), ...
                                          'upperBound', ref_model.ub(refIndex), ...
                                          'objectiveCoef', ref_model.c(refIndex), ...
                                          'reactionName', ref_model.rxnNames{refIndex}, ...
                                          'subSystem', ref_model.subSystems{refIndex}, ...
                                          'geneRule', ref_model.grRules{refIndex});

                    fprintf('Successfully added reaction %s from ref_model to model_1.\n', rxn_add_id);
                else
                    fprintf('Reaction %s not found in ref_model. Addition aborted.\n', rxn_add_id);
                end
            else
                fprintf('Reaction %s already exists in model_1. No addition necessary.\n', rxn_add_id);
            end
        end
    end




    % Option 3: Check for specific reactions
    if ismember(3, option)
        fprintf('\nChecking specified reactions:\n');
        for i = 1:length(checkReactions)
            reaction = checkReactions{i};
            if ismember(reaction, model_1.rxns)
                fprintf('Reaction %s is present in the model.\n', reaction);
            else
                fprintf('Reaction %s is NOT present in the model.\n', reaction);
            end
        end
    end
end