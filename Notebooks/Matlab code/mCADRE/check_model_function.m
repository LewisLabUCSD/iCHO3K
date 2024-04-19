function [requiredMetsStatus, time] = check_model_function(model, varargin)

% Inputs
% - model:
% - optional inputs:
%       'requiredMets', metList
%       'biomass', biomassRxn
%       'media', mediaDef

%% Parse input parameters

if numel(varargin)
    if rem(numel(varargin), 2) == 0
        options = {'requiredMets', 'biomass', 'media'};
        for i = 1:2:numel(varargin)
            argname = varargin{i}; argval = varargin{i + 1};
            option = find(strncmpi(argname, options, numel(argname)));
            if ~isempty(option)
                switch option
                    case 1
                        metList = argval;
                    case 2
                        error('%s option not yet implemented.', argname);
                        % biomassRxn = argval;
                    case 3
                        error('%s option not yet implemented.', argname);
                        % mediaDef = argval;
                end
            else error('Unknown option %s.', argname)
            end
        end
    else error('Incorrect number of input arguments to function %s.', ...
            mfilename);
    end
end

%%
t0 = clock;

% Identify exchange reactions in the model
exRxns = find_ex_rxns(model);

% Turn off uptake of organic metabolites
if exist('mediaDef', 'var')
    model = set_media_ex_bounds(model); % not implemented in this version
else
    model = set_organic_met_bounds(model, exRxns);
end

% Allow uptake of glucose and CO2
model = changeRxnBounds(model, 'EX_glc(e)', -5, 'l');
model = changeRxnBounds(model, 'EX_co2(e)', -1000, 'l');

% Add demand reactions for required metabolites
if exist('metList', 'var')
    [model, requiredRxns] = specify_required_rxns(model, metList);
else requiredRxns = {};
end

if exist('biomassRxn', 'var')
    requiredRxns = [requiredRxns, biomassRxn]; % not implemented in this version
end

inactiveRequired = check_rxn_flux(model, requiredRxns);

requiredMetsStatus = ~numel(inactiveRequired);
time = etime(clock, t0);
