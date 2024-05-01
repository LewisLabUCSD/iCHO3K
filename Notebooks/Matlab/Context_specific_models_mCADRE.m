% Initialize COBRA
initCobraToolbox;

% Path
path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';


% Change COBRA solver to CPLEX or quadMinos
changeCobraSolver('ibm_cplex', 'all');

% Define method and salvageCheck
% method = 1; % fastFVA
method = 2; % fastcc
salvageCheck = 1;

% Reactions names that are manually added to the core reaction
protectedRxns = {'biomass_cho', 'DNAsyn', 'LipidSyn', 'PROTsyn', 'RNAsyn'};

% Load Confidence Scores
confidenceScores = readmatrix(fullfile(path, 'Data/Context_specific_models', 'confidence_scores.csv'));

% Load Model
model = readCbModel(fullfile(path, 'Notebooks', 'iCHO3644_unblocked.mat'));

% Load Ubiquity Scores
load(fullfile(path, 'Data/Context_specific_models', 'UbiquityScores.mat'));
ubiScoreVector = ubiScore(:, 1);

% Define checkFunctionality, eta, and tol
checkFunctionality = 1; % Boolean variable that determines if the model should be able to produce the metabolites associated with the protectedRxns
eta = 1/3; % Trade-off between removing core and zero-expression reactions (default value: 1/3)
tol = 1e-8; % Minimum flux threshold for "expressed" reactions (default 1e-8)

% Execute mCADRE algorithm
[tissueModel, coreRxn, nonCoreRxn, zeroExpRxns, pruneTime, cRes] = mCADRE(model, ubiScoreVector, confidenceScores, protectedRxns, checkFunctionality, eta, tol);
