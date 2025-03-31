% Step 1: Initialize the COBRA toolbox
initCobraToolbox(0);

path = 'C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells';
addpath('C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells\Notebooks\Matlab\Model_Extraction');

model_files = dir('C:\Users\user\Documents\DC\Manual curation_iCHO\Whole-Cell-Network-Reconstruction-for-CHO-cells_origin\Whole-Cell-Network-Reconstruction-for-CHO-cells\Data\Context_specific_models\Reduced models/reduced_model_*.mat');

% Preallocate results
num_models = length(model_files);
dFBA_results(num_models) = struct('time', [], 'profiles', [], 'flux_growth', [], 'flux_prod', [], 'phase_transition', [], 'notes', [], 'condition', []);

% Define common variables
time_range = [0, 10]; % Example time range, adjust as needed
fix_flx = ...; % Define fixed fluxes
fix_comp = ...; % Define fixed components
notes = 'Simulation notes'; % Add any notes if necessary

% Start parallel for loop
for i = 1:num_models
    % Load the current model
    model_file = fullfile(model_files(i).folder, model_files(i).name);
    model = load(model_file);
    
    % Load or define dFBA_data for the current model
    dFBA_data = ...; % Load or define the dFBA data structure
    
    % Run the COSMIC_dFBA simulation
    dFBA_results(i) = COSMIC_dFBA(model, time_range, dFBA_data, fix_flx, fix_comp, notes);
end

% Save results
save('dFBA_results.mat', 'dFBA_results');

% Shut down parallel pool
delete(gcp('nocreate'));