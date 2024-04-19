function E_X = calc_expr_evidence(model, GPRrxns, U_GPR, is_C_H)

E_X = zeros(size(model.rxns));
U_GPR_min = min(U_GPR, [], 2);
for i = 1:numel(model.rxns)
    rxn_GPRs = strmatch(model.rxns{i}, GPRrxns, 'exact');
    if numel(rxn_GPRs)
        E_X(i) = max(U_GPR_min(rxn_GPRs));
    end
end
% For reactions with no corresponding probe in expression data
E_X(isnan(E_X)) = 0;
E_X(is_C_H) = 1;