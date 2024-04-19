function deadEnd_C = check_core_deadends(model, C)

deadEnd_C = [];

for metNo = 1:length(model.mets)
    de = 0;
    
    bothRxns = find(model.S(metNo, :) ~= 0 & model.rev(:)');
    prodRxns = union(find(model.S(metNo,:) > 0), bothRxns);
    consRxns = union(find(model.S(metNo,:) < 0), bothRxns);
    
    % Check for produced-only metabolites
    if (isempty(consRxns))
        de = 1;
    
    % Check for consumed-only metabolites
    elseif (isempty(prodRxns))
        de = 1;
    
    % Check for metabolites both consumed and produced, but only in a single
    % reversible reaction
    else
        if (length(prodRxns) == 1 && ...
                length(consRxns) == 1 && ...
                prodRxns == consRxns)
            de = 1;
        end
    end
    
    % If dead end found, check for overlap with core reactions
    if de == 1
        contains_deadEnd = sum(full(model.S(metNo, :) ~= 0), 1) > 0;
        if any(ismember(model.rxns(contains_deadEnd), C))
            deadEnd_C = model.rxns(contains_deadEnd);
            break;
        end
    end
end

