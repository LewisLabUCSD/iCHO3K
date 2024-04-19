
% Script to prepare main input data files for mcadre models

% Human model: Recon 1 (from BiGG database, May 2011)
clear

load('HR1_CbModel');
load('confidenceScores');

model = format_human_model(model);

% create consistent version for testing
inactiveRxns = check_model_consistency(model);
model_consistent = removeRxns(model, inactiveRxns);

save('../data/humanModel', 'model', 'confidenceScores', 'model_consistent');


% Mouse model: update of iMM1415 (from A. Heinken / I. Thiele, March 2012)
clear

load('MouseModel_New');

model = format_mouse_model(MouseModel_New);
confidenceScores = model.confidenceScores;

% create consistent version for testing
inactiveRxns = check_model_consistency(model);
model_consistent = removeRxns(model, inactiveRxns);

save('../data/mouseModel_nonFunctional', 'model', 'confidenceScores', ...
    'model_consistent');


% Update mouse model to enable generic functionality
model = update_mouse_model(model);
confidenceScores = model.confidenceScores;

% create consistent version for testing
inactiveRxns = check_model_consistency(model);
model_consistent = removeRxns(model, inactiveRxns);

save('../data/mouseModel', 'model', 'confidenceScores', 'model_consistent');
