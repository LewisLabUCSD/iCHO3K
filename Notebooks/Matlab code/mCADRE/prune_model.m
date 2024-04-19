function [PM, cRes] = ...
    prune_model(GM, P, C, Z, eta, precursorMets, salvageCheck, method, cutoff)

% Initialize variables
R_G = GM.rxns;
PM = GM;
R_P = R_G;

NC_removed = 0; C_removed = 0;
cRes = zeros(3000, 1);
count = 1;

if nargin < 9
    cutoff = Inf;
end

if nargin < 8
    method = 1; % fastFVA
end

if nargin < 7
    salvageCheck = 1; % assume non-hepatic tissue, where nucleotides are not
                      % expected to be produced de novo
end

while numel(P) && count < cutoff % for testing
    display(['Reaction no. ', num2str(count)])
    r = P(1);
    display(['Attempting to remove reaction ', r{:}, '...'])
    modelR = removeRxns(PM, r);

    % First check precursor production; if this test fails, no need to
    % check model consistency with FVA (time-saving step)
    rStatus = check_model_function(modelR, ...
        'requiredMets', precursorMets);
    
    % If specified, check the salvage pathway as well
    if salvageCheck
        rSalvage = check_salvage_path(modelR);
        rStatus = rStatus && rSalvage;
    end
    
    if rStatus

        % Check for inactive reactions after removal of r
        inactive_G = check_model_consistency(PM, method, r);

        inactive_C = intersect(inactive_G, C);
        inactive_NC = setdiff(inactive_G, inactive_C);

        % Remove reactions with zero expression (previously penalized in
        % rank_reactions) and corresponding inactive core reactions, only if
        % sufficiently more non-core reactions are removed
        if ismember(r, Z)
            display('Zero-expression evidence for reaction...')
            
            % Check model function with all inactive reactions removed
            modelTmp = removeRxns(PM, inactive_G);
            tmpStatus = check_model_function(modelTmp, ...
                'requiredMets', precursorMets);
            
            % If specified, check the salvage pathway as well
            if salvageCheck
                tmpSalvage = check_salvage_path(modelTmp);
                tmpStatus = tmpStatus && tmpSalvage;
            end

            if (numel(inactive_C) / numel(inactive_NC) <= eta) && tmpStatus
                R_P = setdiff(R_P, inactive_G);
                PM = removeRxns(PM, inactive_G);
                P(ismember(P, inactive_G)) = [];
                NC_removed = NC_removed + numel(inactive_NC);
                C_removed = C_removed + numel(inactive_C);
                num_removed = NC_removed + C_removed;
                display('Removed all inactive reactions')

                % result = -1.x indicates that reaction r had zero
                % expression evidence and was removed along with any
                % consequently inactivated reactions; x indicates the number of
                % core reactions removed
                if numel(inactive_C) > 100
                    removed_C_indicator = numel(inactive_C) / 100;
                else removed_C_indicator = numel(inactive_C) / 10;
                end
                result = -1 - removed_C_indicator;
            else
                % Note: no reactions (core or otherwise) are actually
                % removed in this step, but it is necessary to update the
                % total number of removed reactions to avoid errors below
                num_removed = NC_removed + C_removed;
                P(1) = [];
                display('No reactions removed')

                % result = 1.x indicates that no reactions were removed
                % because removal of r either led to a ratio of inactivated
                % core vs. non-core reactions above the specified threshold
                % eta (x = 1) or the removal of r and consequently
                % inactivated reactions prevented production of required
                % metabolites (x = 0)
                result = 1 + tmpStatus / 10;
            end

        % If reaction has expression evidence, only attempt to remove
        % inactive non-core reactions
        else
            % Check model function with non-core inactive reactions removed
            modelTmp = removeRxns(PM, inactive_NC);
            tmpStatus = check_model_function(modelTmp, ...
                'requiredMets', precursorMets);
            
            % If specified, check the salvage pathway as well
            if salvageCheck
                tmpSalvage = check_salvage_path(modelTmp);
                tmpStatus = tmpStatus && tmpSalvage;
            end

            if numel(inactive_C) == 0 && tmpStatus
                R_P = setdiff(R_P, inactive_NC);
                PM = removeRxns(PM, inactive_NC);
                P(ismember(P, inactive_NC)) = [];
                NC_removed = NC_removed + numel(inactive_NC);
                num_removed = NC_removed + C_removed;
                display('Removed non-core inactive reactions')

                % result = -2 indicates that reaction r had expression.
                % evidence and was removed along with (only) non-core
                % inactivated reactions; x indicates the number of
                % core reactions removed (should be zero!)
                if numel(inactive_C) > 100
                    removed_C_indicator = numel(inactive_C) / 100;
                else removed_C_indicator = numel(inactive_C) / 10;
                end
                result = -2 - removed_C_indicator;
            else
                num_removed = NC_removed + C_removed;
                P(1) = [];
                display('No reactions removed')

                % result = 2.x indicates that no reactions were removed
                % because removal of r either led to inactivated core
                % reactions (x = 1) or the removal of r and consequently
                % inactivated reactions prevented production of required
                % metabolites (x = 0)
                result = 2 + tmpStatus / 10;
            end
        end
    else
        num_removed = NC_removed + C_removed;
        P(1) = [];

        % result = 3 indicates that no reactions were removed because
        % removal of r by itself prevented production of required
        % metabolites
        result = 3;
    end
    
    cRes(count) = result;
    count = count + 1;
    display(sprintf(['Num. removed: ', num2str(num_removed), ...
        ' (', num2str(C_removed), ' core, ', ...
        num2str(NC_removed), ' non-core); ', ...
        'Num. remaining: ', num2str(numel(P)), '\n']))
end
cRes(count:end) = [];
