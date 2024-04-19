function organicExRxns = find_organic_ex_rxns(model, exRxns)

% Note: should add a warning if the metFormulas field is empty

% Organic metabolites are defined as those containing carbon (C) and hydrogen
% (H); these are identified by checking molecular formulas
carbonMets = ~cellfun('isempty', regexp(model.metFormulas, 'C'));
hydrogenMets = ~cellfun('isempty', regexp(model.metFormulas, 'H'));
is_organic = carbonMets & hydrogenMets;
organicMets = model.mets(is_organic);

organicRxns = findRxnsFromMets(model, organicMets);
organicExRxns = intersect(organicRxns, exRxns);

% The following reactions exchange organic metabolites (e.g., R-groups that
% comprise lipid tails), but don't contain H in their specified formulas; OR,
% as in the case of Tyr-ggn, include protein compounds
organicExRxns = [organicExRxns; ...
    'EX_Rtotal(e)'; 'EX_Rtotal2(e)'; 'EX_Rtotal3(e)'; 'EX_Tyr_ggn(e)']; % ; ...
    % 'UP_Tyr_ggn[c]']; % This rxn is not in Recon 1 - may be something specific
    % to Recon 2...