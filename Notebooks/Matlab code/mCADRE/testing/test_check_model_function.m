
load('humanModel');
load('precursorMets');
changeCobraSolver('glpk');

%% Test find_ex_rxns

display('## Testing find_ex_rxns...');
try
    exRxns = find_ex_rxns(model);
    display(['PASS...', ...
        'Function find_ex_rxns ran without error']);
catch err
    display(['FAIL...', ...
        'Function find_ex_rxns was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

%% Test set_media_ex_bounds

% mediaDef option not implemented in this version


%% Test set_organic_met_bounds

check = false(2, 1);

% Test find_organic_ex_rxns
display('## Testing find_organic_ex_rxns...');

try
    organicExRxns = find_organic_ex_rxns(model, exRxns);
    check(1) = true;
    display(['PASS...', ...
        'Function find_ex_rxns ran without error']);
catch err
    display(['FAIL...', ...
        'Function find_ex_rxns was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Test set_organic_met_bounds
display('## Testing set_organic_met_bounds...');

try
    testModel = set_organic_met_bounds(model, exRxns);
    check(2) = true;
    display(['PASS...', ...
        'Function set_organic_met_bounds ran without error']);
catch err
    display(['FAIL...', ...
        'Function set_organic_met_bounds was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check that all organicExRxn bounds were set correctly
display('## Testing output of set_organic_met_bounds...');

if all(check)
    is_organic_ex = ismember(testModel.rxns, organicExRxns);
    if all(testModel.lb(is_organic_ex) == 0)
        display(['PASS...', ...
            'Exchange reaction bounds were set correctly']);
    else
        display(['FAIL...', ...
            'At least one organic exchange reaction lower-bound was ', ...
            'not correctly set to 0.']);
    end
else
    display(['FAIL...', ...
        'At least one of the above tests did not pass.'])
end
display('---');
    

%% Test specify_required_rxns

display('## Testing specify_required_rxns');

metList = precursorMets;
try
    [testModel, requiredRxns] = specify_required_rxns(model, metList);
    display(['PASS...', ...
        'Function specify_required_rxns ran without error']);
catch err
    display(['FAIL...', ...
        'Function specify_required_rxns was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');


%% Test check_rxn_flux

display('## Checking for success on preceding tests...');

try
    metList = precursorMets;
    exRxns = find_ex_rxns(model);
    testModel = set_organic_met_bounds(model, exRxns);
    [testModel, requiredRxns] = specify_required_rxns(testModel, metList);
    
    % Allow uptake of glucose and CO2
    testModel = changeRxnBounds(testModel, 'EX_glc(e)', -5, 'l');
    testModel = changeRxnBounds(testModel, 'EX_co2(e)', -1000, 'l');
catch err
    display(['FAIL...', ...
        'Cannot test check_rxn_flux because one of the preceding functions', ...
        'returned the error:']);
    display(['> ', err.message]);
end
display('---');

% Test whether Recon 1 passes functional check out of the box (it should)
display('## Testing check_rxn_flux with functional input...');

try
    inactiveRequired = check_rxn_flux(testModel, requiredRxns);
    display(['PASS...', ...
        'Function check_rxn_flux ran without error']);
    if ~numel(inactiveRequired)
        display(['PASS...', ...
            'Check for functionality returns expected result']);
    else
        display(['FAIL...', ...
            'Check for functionality returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_rxn_flux was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Test if Recon 1 passes check when you 'break' it (turn off glucose)
display('## Testing check_rxn_flux with glucose uptake off...');

try
    testModel_noGlc = changeRxnBounds(testModel, 'EX_glc(e)', 0, 'l');
    inactiveRequired = check_rxn_flux(testModel_noGlc, requiredRxns);
    display(['PASS...', ...
        'Function check_rxn_flux ran without error']);
    if numel(inactiveRequired)
        display(['PASS...', ...
            'Check for functionality returns expected result']);
    else
        display(['FAIL...', ...
            'Check for functionality returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_rxn_flux was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Test if Recon 1 passes check when you 'break' it (output blocked)
display('## Testing check_rxn_flux with one output blocked...');

try
    testModel_block = changeRxnBounds(testModel, 'DM_3pg[c]', 0, 'u');
    inactiveRequired = check_rxn_flux(testModel_block, requiredRxns);
    display(['PASS...', ...
        'Function check_rxn_flux ran without error']);
    if numel(inactiveRequired)
        display(['PASS...', ...
            'Check for functionality returns expected result']);
    else
        display(['FAIL...', ...
            'Check for functionality returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_rxn_flux was terminated with the error:']);
    display(['> ', err.message]);
end  
display('---');


%% Test check_model_function

% Test with default inputs
display('## Testing check_model_function with default inputs...')
try
    genericStatus = check_model_function(model, ...
        'requiredMets', precursorMets);
    display(['PASS...', ...
        'Function check_model_function ran without error']);
    
    if genericStatus
        display(['PASS...', ...
            'Check for functionality returns expected result']);
    else
        display(['FAIL...', ...
            'Check for functionality returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_model_function was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Test with reactions removed that will prevent function
display('## Testing check_model_function with glucose reactions removed...')

glcRxns = findRxnsFromMets(model, 'glc-DASH-D[e]');
glcRxns = setdiff(glcRxns, 'EX_glc(e)');
modelR = removeRxns(model, glcRxns);

try
    genericStatus = check_model_function(modelR, ...
        'requiredMets', precursorMets);
    display(['PASS...', ...
        'Function check_model_function ran without error']);
    
    if ~genericStatus
        display(['PASS...', ...
            'Check for functionality returns expected result']);
    else
        display(['FAIL...', ...
            'Check for functionality returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_model_function was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Test with 'biomass' flag
display('## Testing check_model function with non-implemented option...');

try
    genericStatus = check_model_function(model, ...
        'requiredMets', precursorMets, ...
        'biomass', 'dummy_rxn');
    display(['PASS...', ...
        'Function check_model_function ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_model_function was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');   

% Test with 'biomass' flag
display('## Testing check_model function with non-implemented option...');

try
    genericStatus = check_model_function(model, ...
        'requiredMets', precursorMets, ...
        'media', 'dummy_media');
    display(['PASS...', ...
        'Function check_model_function ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_model_function was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');


%% Test check_salvage_path

% Note: not directly part of check_model_function, but closely related

% Test with default inputs
display('## Testing check_salvage_path with default inputs...')
try
    salvageStatus = check_salvage_path(model);
    display(['PASS...', ...
        'Function check_salvage_path ran without error']);
    
    if salvageStatus
        display(['PASS...', ...
            'Check for salvage pathway returns expected result']);
    else
        display(['FAIL...', ...
            'Check for salvage pathway returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_salvage_path was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');


%% Test check_model_function with mouse

load('mouseModel');
load('precursorMets');
changeCobraSolver('glpk');

% Test with default inputs
display('## Testing check_model_function (mouse) with default inputs...')
try
    genericStatus = check_model_function(model, ...
        'requiredMets', precursorMets);
    display(['PASS...', ...
        'Function check_model_function ran without error']);
    
    if genericStatus
        display(['PASS...', ...
            'Check for functionality returns expected result']);
    else
        display(['FAIL...', ...
            'Check for functionality returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_model_function was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Test with reactions removed that will prevent function
display('## Testing check_model_function with glucose reactions removed...')

glcRxns = findRxnsFromMets(model, 'glc_D[e]');
glcRxns = setdiff(glcRxns, 'EX_glc(e)');
modelR = removeRxns(model, glcRxns);

try
    genericStatus = check_model_function(modelR, ...
        'requiredMets', precursorMets);
    display(['PASS...', ...
        'Function check_model_function ran without error']);
    
    if ~genericStatus
        display(['PASS...', ...
            'Check for functionality returns expected result']);
    else
        display(['FAIL...', ...
            'Check for functionality returns unexpected result']);
    end
catch err
    display(['FAIL...', ...
        'Function check_model_function was terminated with the error:']);
    display(['> ', err.message]);
end