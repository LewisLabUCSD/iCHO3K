function model = update_mouse_model(model)

% The model contained in `MouseModel_New.mat` is an updated version of iMM1415,
% published by Sirudsson et al. in 2010 [1]. This model, when tested for
% functionality with `check_model_function()`, fails to produce some required
% precursor metabolites. In her thesis work, Chunjing Wang identified two
% reactions (with genetic evidence) that enabled the production of blocked
% non-essential amino acids [2]. These reactions have been added here, such that
% the `mouseModel` passes the check for generic functionality, prior to removal
% of any reactions.

% 1. M Sigurdsson, et al. A detailed genome-wide reconstruction of mouse
%    metabolism based on human Recon 1. BMC Syst Biol (2010)
% 2. C Wang. Characterizing and analyzing disease-related omics data using
%    network modeling approaches. Dissertation, University of Illinois (2014)

% Add glutamate synthesis reactions
[~, model] = evalc(['addReaction(model, ''GLUDxm'',', ...
    '{''h2o[m]''; ''nad[m]''; ''glu_L[m]''; ', ...
    '''h[m]''; ''akg[m]''; ''nadh[m]''; ''nh4[m]''},', ...
    '[-1; -1; -1; 1; 1; 1; 1], 1, -10000, 10000, 0', ...
    ''''', ''14661.1 or 104277.1'');']);

[~, model] = evalc(['addReaction(model, ''GLUDym'',', ...
    '{''h2o[m]''; ''nadp[m]''; ''glu_L[m]''; ', ...
    '''h[m]''; ''akg[m]''; ''nadph[m]''; ''nh4[m]''},', ...
    '[-1; -1; -1; 1; 1; 1; 1], 1, -10000, 10000, 0', ...
    ''''', ''14661.1 or 104277.1'');']);
