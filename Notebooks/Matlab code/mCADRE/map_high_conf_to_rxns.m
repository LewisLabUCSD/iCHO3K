function is_C_H = map_high_conf_to_rxns(model, GPRmat, GPRrxns, C_H_genes)
    if isempty(C_H_genes)
        is_C_H = false(size(model.rxns));  % No high-confidence genes provided, all reactions are not high-confidence
        return;
    end

    % Ensure all elements in C_H_genes are strings (this should be redundant now)
    if isnumeric(C_H_genes) || ischar(C_H_genes)
        C_H_genes = cellstr(num2str(C_H_genes));
    end

    % Mapping logic
    C_H_GPR = double(ismember(GPRmat, C_H_genes));
    C_H_GPR(cellfun('isempty', GPRmat)) = NaN;  % Handle empty gene-reaction associations
    C_H_GPR_min = min(C_H_GPR, [], 2);  % Find minimum across rows
    is_C_H = false(size(model.rxns));
    for i = 1:numel(model.rxns)
        rxn_GPRs = find(strcmp(model.rxns{i}, GPRrxns));
        if ~isempty(rxn_GPRs)
            is_C_H(i) = any(C_H_GPR_min(rxn_GPRs) == 1);
        end
    end
    is_C_H = logical(is_C_H);
end

