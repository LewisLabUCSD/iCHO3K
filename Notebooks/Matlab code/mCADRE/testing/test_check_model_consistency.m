
load('humanModel');
changeCobraSolver('glpk');

%% Test check_core_deadends

% Note: this function isn't actually needed for mCADRE, but was originally
% included to help speed up my code for running MBA. As mCADRE allows for the
% removal of core reactions under some scenarios, this function shouldn't be
% called with check_model_consistency (i.e., deCheck = 0).

% Test with random core
display('## Testing check_core_deadends with random core...');

rng(0) % set the seed to enable reproduciblity
rand_idx = randperm(length(model.rxns), 500);
C = model.rxns(rand_idx);

try
    deadEnd_C = check_core_deadends(model, C);
    display(['PASS...', ...
        'Function check_core_deadends ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_core_deadends was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 1 dead-end-containing core reaction
display('## Checking output of check_core_deadends with random core...');

if exist('deadEnd_C', 'var')
    if numel(deadEnd_C) == 1
        display(['PASS...', ...
            'Function check_core_deadend returns the expected result']);
    else
        display(['FAIL...', ...
            'Function check_core_deadend returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function check_core_deadend encountered error, cannot check output'])
end

% Test with small core which should contain no dead-end metabolites
display('## Testing check_core_deadends with small core...');

C = {'EX_glc(e)'};

try
    deadEnd_C = check_core_deadends(model, C);
    display(['PASS...', ...
        'Function check_core_deadends ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_core_deadends was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find no dead-end-containing core reactions
display('## Checking output of check_core_deadends with small core...');

if exist('deadEnd_C', 'var')
    if numel(deadEnd_C) == 0
        display(['PASS...', ...
            'Function check_core_deadend returns the expected result']);
    else
        display(['FAIL...', ...
            'Function check_core_deadend returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function check_core_deadend encountered error, cannot check output'])
end


%% Test find_inactive_rxns

% Test with fastFVA
display('## Testing find_inactive_rxns with fastFVA...');

try
    inactiveRxns = find_inactive_rxns(model, 1);
    display(['PASS...', ...
        'Function find_inactive_rxns ran without error']);
catch err
    display(['FAIL...', ...
        'Function find_inactive_rxns was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 1273 inactive reactions in Recon 1
display('## Checking output of find_inactive_rxns with fastFVA...');

if exist('inactiveRxns', 'var')
    if numel(inactiveRxns) == 1273
        display(['PASS...', ...
            'Function find_inactive_rxns returns the expected result']);
    else
        display(['FAIL...', ...
            'Function find_inactive_rxns returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function find_inactive_rxns encountered error, cannot check output'])
end

% Test with fastcc
display('## Testing find_inactive_rxns with fastcc...');

clear('inactiveRxns')
try
    inactiveRxns = find_inactive_rxns(model, 2);
    display(['PASS...', ...
        'Function find_inactive_rxns ran without error']);
catch err
    display(['FAIL...', ...
        'Function find_inactive_rxns was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 1273 inactive reactions in Recon 1
display('## Checking output of find_inactive_rxns with fastFVA...');

if exist('inactiveRxns', 'var')
    if numel(inactiveRxns) == 1273
        display(['PASS...', ...
            'Function find_inactive_rxns returns the expected result']);
    else
        display(['FAIL...', ...
            'Function find_inactive_rxns returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function find_inactive_rxns encountered error, cannot check output'])
end


%% Test check_model_consistency

% Test with fastFVA, default parameters
method = 1;
r = [];
deCheck = 0;
C = {};

display('## Testing check_model_consistency with default inputs...');

try
    [inactiveRxns, time, result] = check_model_consistency(model, ...
        method, r, deCheck, C);
    display(['PASS...', ...
        'Function check_model_consistency ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_model_consistency was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 1273 inactive reactions in Recon 1
display('## Checking output of check_model_consistency with default inputs...');

if exist('inactiveRxns', 'var')
    if numel(inactiveRxns) == 1273
        display(['PASS...', ...
            'Function check_model_consistency returns the expected result']);
    else
        display(['FAIL...', ...
            'Function check_model_consistency returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function check_model_consistency encountered error, no output'])
end

% Test with fastFVA, consistent model
method = 1;
r = [];
deCheck = 0;
C = {};

display('## Testing check_model_consistency with consistent model...');

try
    [inactiveRxns, time, result] = check_model_consistency(model_consistent, ...
        method, r, deCheck, C);
    display(['PASS...', ...
        'Function check_model_consistency ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_model_consistency was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 0 inactive reactions in Recon 1
display('## Checking output of check_model_consistency with default inputs...');

if exist('inactiveRxns', 'var')
    if numel(inactiveRxns) == 0
        display(['PASS...', ...
            'Function check_model_consistency returns the expected result']);
    else
        display(['FAIL...', ...
            'Function check_model_consistency returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function check_model_consistency encountered error, no output'])
end

% Test with fastFVA, using consistent model + removed reaction
r = 'ACCOAC';

display('## Testing check_model_consistency with reaction removed');

try
    [inactiveRxns, time, result] = check_model_consistency(model_consistent, ...
        method, r, deCheck, C);
    display(['PASS...', ...
        'Function check_model_consistency ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_model_consistency was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 13 inactive reactions in Recon 1
display('## Checking output of check_model_consistency, reaction removed...');

if exist('inactiveRxns', 'var')
    if numel(inactiveRxns) == 13
        display(['PASS...', ...
            'Function check_model_consistency returns the expected result']);
    else
        display(['FAIL...', ...
            'Function check_model_consistency returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function check_model_consistency encountered error, no output'])
end


%%
% Test with fastcc
display('## Testing check_model_consistency with fastcc...');

clear('inactiveRxns')
try
    [inactiveRxns, time, result] = check_model_consistency(model, 2);
    display(['PASS...', ...
        'Function check_model_consistency ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_model_consistency was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 1273 inactive reactions in Recon 1
display('## Checking output of check_model_consistency with fastcc...');

if exist('inactiveRxns', 'var')
    if numel(inactiveRxns) == 1273
        display(['PASS...', ...
            'Function check_model_consistency returns the expected result']);
    else
        display(['FAIL...', ...
            'Function check_model_consistency returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function check_model_consistency encountered error, no output'])
end

% Test with fastcc, consistent model
method = 2;
r = [];
deCheck = 0;
C = {};

display('## Testing check_model_consistency with fastcc, consistent model...');

try
    [inactiveRxns, time, result] = check_model_consistency(model_consistent, ...
        method, r, deCheck, C);
    display(['PASS...', ...
        'Function check_model_consistency ran without error']);
catch err
    display(['FAIL...', ...
        'Function check_model_consistency was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs; should find 0 inactive reactions in Recon 1
display('## Checking output of check_model_consistency with fastcc...');

if exist('inactiveRxns', 'var')
    if numel(inactiveRxns) == 0
        display(['PASS...', ...
            'Function check_model_consistency returns the expected result']);
    else
        display(['FAIL...', ...
            'Function check_model_consistency returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function check_model_consistency encountered error, no output'])
end