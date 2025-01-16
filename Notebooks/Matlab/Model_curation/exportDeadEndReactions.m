function exportDeadEndReactions(newModel,oldModel, addedRxns, useBiGGIDorNames, useTimestampedFilename)
%%
% This function exports information about reactions added to a metabolic model
% and identifies dead-end metabolites in the updated model. It can optionally use
% BiGG IDs or reaction names instead of default reaction IDs and save the output
% with a timestamp.
%
% Inputs:
%   - newModel: The updated metabolic model after gap filling (structure).
%   - oldModel: The original metabolic model before gap filling (structure).
%   - addedRxns: Cell array of IDs of reactions added to the newModel.
%
% Optional Inputs:
%   - useBiGGIDorNames: Logical (true/false). If true, uses BiGG IDs (if possible) or reaction names
%     instead of default reaction IDs. Default is false.
%   - useTimestampedFilename: Logical (true/false). If true, the exported Excel file
%     will include a timestamp in its filename. Default is false.
%
% Outputs:
%   - An Excel file named 'dead_ends_reactions.xlsx' or with a timestamp(optional), containing
%     two sheets: 'Added Reactions' and 'Filtered Reactions'. These sheets document
%     reactions added to the model and information about dead-end metabolites and their
%     related reactions, respectively.
% Proceeding code example
% % Initialize the COBRA Toolbox without updating the path
% initCobraToolbox(false);
% 
% % Import the original metabolic model
% oldModel = importModel('~filepath\model_name.xml');
% 
% % Import additional models for gap filling
% Recon_model_raven = importModel('Recon3D.xml');
% iCHO_model_raven = importModel('iCHO_v1.xml');
% iMM_model_raven = importModel('iMM1415.xml');
% 
% % Combine additional models into a cell array
% RCM_model_raven = {Recon_model_raven, iCHO_model_raven, iMM_model_raven};
% 
% % Perform gap filling on the original model using the additional models
% [RCM_newConnected, RCM_cannotConnect, RCM_addedRxns, newModel, RCM_exitFlag] = fillGaps(oldModel, RCM_model_raven);
% 
% % Now, you can use the `exportDeadEndReactions` function with the prepared inputs
% % Example:
% % exportDeadEndReactions(newModel, oldModel, RCM_addedRxns, false, true);
%%
% Default values for optional parameters
    if nargin < 4
        useBiGGIDorNames = false;
    end
    if nargin < 5
        useTimestampedFilename = false;
    end

    % Print added reaction formulas
    added_reactionFormulas = cell(length(addedRxns), 2);
    for i_len = 1:length(addedRxns)
        formula_tmp = printRxnFormula(newModel, {addedRxns{i_len}}, false, true, false, 1);
        reactionID = {addedRxns{i_len}};
        
        if useBiGGIDorNames
            if isfield(newModel, 'rxnBiGGID') && ~isempty(newModel.rxnBiGGID{i_len})
                reactionID = {newModel.rxnBiGGID{i_len}};
            elseif isfield(newModel, 'rxnNames') && ~isempty(newModel.rxnNames{i_len})
                reactionID = {newModel.rxnNames{i_len}};
            end
        end
        
        added_reactionFormulas{i_len, 1} = reactionID{1, 1};
        added_reactionFormulas{i_len, 2} = formula_tmp{1}; % Store the formula
    end

    addedFormulasTable = cell2table(added_reactionFormulas, 'VariableNames', {'RxnID', 'Formula'});

    % Find existing dead-end metabolites
    DeadEnds = newModel.mets(detectDeadEnds(oldModel, false));
 
    % Initialize variables for storing filtered reactions
    filteredFormulas = {'DeadEndMetaboliteID', 'RelatedReactions'};
    
    for i = 1:length(DeadEnds)
        met = DeadEnds{i};
        len = length(met);
        selectedDeadEnd = met(1:len-1); % Remove last character
        
        % Initialize temporary array for storing related reactions
        tempReactions = {};
        for j = 1:height(addedFormulasTable)
            formula = addedFormulasTable{j, 2}{1}; % Extract reaction formula
            % Check if current dead-end metabolite is included in the reaction formula
            if contains(formula, selectedDeadEnd)
                % Include rxnID and reaction formula in tempReactions
                rxnID = addedFormulasTable{j, 1}{1};
                tempReactions{end+1} = sprintf('%s:%s', rxnID, formula);
            end
        end
        
        % Add to filteredFormulas if tempReactions is not empty
        if ~isempty(tempReactions)
            relatedReactionsStr = strjoin(tempReactions, '; ');
            filteredFormulas = [filteredFormulas; {DeadEnds{i}, relatedReactionsStr}];
        end
    end

    % Convert filtered reactions to a table
    filteredFormulasTable = cell2table(filteredFormulas(2:end,:), 'VariableNames', filteredFormulas(1,:));

    % Define the Excel file name
    if useTimestampedFilename
        timestamp = datestr(now, 'mm-dd-HH-MM');
        filename = sprintf('dead_ends_reactions_%s.xlsx', timestamp);
    else
        filename = 'dead_ends_reactions.xlsx';
    end

    % Save the table to an Excel file
    writetable(addedFormulasTable, filename, 'Sheet', 'Added Reactions');
    writetable(filteredFormulasTable, filename, 'Sheet', 'Filtered Reactions');
end
