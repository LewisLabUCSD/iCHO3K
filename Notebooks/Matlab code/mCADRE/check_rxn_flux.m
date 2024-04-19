function inactiveRequired = check_rxn_flux(model, requiredRxns)

% This function uses the heuristic speed-up proposed by Jerby et al. in the MBA
% paper for performing a pseudo-FVA calculation.

rxnList = requiredRxns;
inactiveRequired = [];
while numel(rxnList)
    numRxnList = numel(rxnList);
    % model.rxns(strmatch('biomass_', model.rxns)); % not implemented
    model = changeObjective(model, rxnList);

    % Maximize all
    FBAsolution = optimizeCbModel(model, 'max');

    optMax = FBAsolution.x;
    % If no solution was achieved when trying to maximize all reactions, skip
    % the subsequent step of checking individual reactions
    if isempty(optMax)
        inactiveRequired = 1;
        break;
    end
    requiredFlux = optMax(ismember(model.rxns, requiredRxns));
    activeRequired = requiredRxns(abs(requiredFlux) >= 1e-6);
    rxnList = setdiff(rxnList, activeRequired);

    numRemoved = numRxnList - numel(rxnList);

    if ~numRemoved
        randInd = randperm(numel(rxnList));
        i = rxnList(randInd(1));
        model = changeObjective(model, i);

        % Maximize reaction i
        FBAsolution = optimizeCbModel(model, 'max');
        optMax = FBAsolution.f;
        if isempty(optMax)
            inactiveRequired = union(inactiveRequired, i);
            break;
        end
        if abs(optMax) < 1e-6
            inactiveRequired = union(inactiveRequired, i);
            break;
        end

        rxnList = setdiff(rxnList, i);
    end
end
