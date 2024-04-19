load('testInputs');
load('humanModel');
load('precursorMets');
changeCobraSolver('glpk');

%% Run rank_reactions

% Not really a better way (that is obvious right now) to test prune_model
try 
    [GM, C, NC, P, Z, model_C] = ... 
         rank_reactions(model, G, U, confidenceScores, C_H_genes);
     display(['READY...', ...
         'Function rank_reactions ran without error; inputs should be', ...
         'for testing prune_model']);
catch err
    display(['FAIL...', ...
        'Function rank_reactions was terminated with the error:']);
    display(['> ', err.message]);
    display('Further testing of rank_reactions is recommended')
end
display('---');

%% Test prune_model

display('## Testing prune_model...');

eta = 1/3;
method = 1;
salvageCheck = 1;

try
    cutoff = 3;
    [~, PM, cRes] = evalc( ...
        ['prune_model(GM, P, C, Z, eta, precursorMets, salvageCheck,', ...
        'method, cutoff);']);
    display(['PASS...', ...
        'Function prune_model ran without error']);
catch err
    display(['FAIL...', ...
        'Function prune_model was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

% Check outputs for zero-expression reactions
display('## Checking outputs with zero-expression reactions...')

try
    cutoff = 11;
    Ptest = P(1:10);
    [~, PM_Z, cRes] = evalc( ...
        ['prune_model(GM, Ptest, C, Z, eta, precursorMets, salvageCheck,', ...
        'method, cutoff);']);
    display(['PASS...', ...
        'Function prune_model ran without error']);
catch err
    display(['FAIL...', ...
        'Function prune_model was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

if exist('PM_Z', 'var')
    resCheck = cRes(cRes < 0);
    if mod(sum(resCheck), 1) == 0.5
        display(['PASS...', ...
            'Function prune_model returns the expected result']);
    else
        display(['FAIL...', ...
            'Function prune_model returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function prune_model encountered error, cannot check output'])
end
display('---');

% Check outputs for expression evidence reactions
display(['## Checking outputs for reactions with missing or positive ', ...
    'expression...']);

try
    cutoff = 11;
    Ptest = P(70:80);
    [~, PM_nZ, cRes] = evalc( ...
        ['prune_model(GM, Ptest, C, Z, eta, precursorMets, salvageCheck,', ...
        'method, cutoff);']);
    display(['PASS...', ...
        'Function prune_model ran without error']);
catch err
    display(['FAIL...', ...
        'Function prune_model was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');

if exist('PM_nZ', 'var')
    resCheck = cRes(cRes < 0);
    if mod(sum(resCheck), 1) == 0
        display(['PASS...', ...
            'Function prune_model returns the expected result']);
    else
        display(['FAIL...', ...
            'Function prune_model returns unexpected result']);
    end
else
    display(['FAIL...', ...
        'Function prune_model encountered error, cannot check output'])
end


%% Test prune_model with fastcc

display('## Testing prune_model with fastcc...');

eta = 1/3;
method = 2;
salvageCheck = 1;

try
    cutoff = 3;
    [~, PM, cRes] = evalc( ...
        ['prune_model(GM, P, C, Z, eta, precursorMets, salvageCheck,', ...
        'method, cutoff);']);
    display(['PASS...', ...
        'Function prune_model ran without error']);
catch err
    display(['FAIL...', ...
        'Function prune_model was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');


%% Test remove_unused_rxns

% Note: not directly part of prune_model, but closely related

% Test with default inputs
display('## Testing remove_unused_rxns with test pruning outputs...')
try
    load('mcadre_results');
    
    PMtest = remove_unused_genes(PM);
    display(['PASS...', ...
        'Function remove_unused_rxns ran without error']);
    
    if exist('PMtest', 'var')
        if numel(PMtest.genes) == 1282
            display(['PASS...', ...
                'Function remove_unused_rxns returns expected result']);
        else
            display(['FAIL...', ...
                'Function remove_unused_rxns returns unexpected result']);
        end
    else
        display(['FAIL...', ...
            'Function remove_unused_rxns returned error, cannot check output']);
    end
catch err
    display(['FAIL...', ...
        'Function remove_unused_rxns was terminated with the error:']);
    display(['> ', err.message]);
end
display('---');