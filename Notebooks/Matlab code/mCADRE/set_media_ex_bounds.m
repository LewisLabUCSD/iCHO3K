function model = set_media_ex_bounds(model, mediaDef)

% Reset all exchanges to zero
model = changeRxnBounds(model,exRxns, 0, 'b');

% Set specified exchanges to specified bounds
[~, modelIdx, mediaIdx] = intersect(model.rxns, mediaDef.rxns);
model.lb(modelIdx) = mediaDef.lb(mediaIdx);
model.ub(modelIdx) = mediaDef.ub(mediaIdx);
