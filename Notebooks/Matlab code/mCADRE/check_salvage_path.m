function [salvageStatus, time] = check_salvage_path(model)

% When the model is allowed to use PRPP and guanine or hypoxanthine, test if it
% can make GMP or IMP. This is the salvage pathway that non-hepatic tissues use
% for purine synthesis. Not useful when the tissue is known to make purines de
% novo.

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

% Add PRPP sink reaction for subsequent checks
[~, model] = evalc('addSinkReactions(model, {''prpp[c]''}, -5, 5);');

% Check production of GMP:

% Allow uptake of guanine
model_gmp = changeRxnBounds(model, 'EX_gua(e)', -5, 'l');

% Add demand reaction for GMP
[~, model_gmp, gmp_dm] = evalc('addDemandReaction(model_gmp, ''gmp[c]'');');
model_gmp = changeObjective(model_gmp, gmp_dm, 1);
sol = optimizeCbModel(model_gmp);
status_gmp = sol.f > 1e-6;

% Check production of IMP:

% Allow uptake of hypoxanthine
model_imp = changeRxnBounds(model, 'EX_hxan(e)', -5, 'l');

% Add demand reaction for IMP
[~, model_imp, imp_dm] = evalc('addDemandReaction(model_imp, ''imp[c]'');');
model_imp = changeObjective(model_imp, imp_dm, 1);
sol = optimizeCbModel(model_imp);
status_imp = sol.f > 1e-6;

salvageStatus = status_gmp && status_imp;

time = etime(clock, t0);