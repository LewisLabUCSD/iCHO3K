function inactiveRxns = find_inactive_rxns(model, method)

% fastFVA is the default method
if nargin < 2
    method = 1;
end

% Check for inactive reactions with either fastFVA or fastcc
if method == 1
    display('Checking all reactions (fastFVA)...')
    model.c(logical(model.c)) = 0;
    [optMin, optMax] = fastFVA(model, 0, 'max', 'glpk');
    is_inactive = (abs(optMax) < 1e-6) & (abs(optMin) < 1e-6);
    inactiveRxns = model.rxns(is_inactive);

else % otherwise, use FASTCC
    display('Checking all reactions (FASTCC)...')
    is_active = fastcc(model, 1e-4);
    inactiveRxns = setdiff(model.rxns, model.rxns(is_active));
end
