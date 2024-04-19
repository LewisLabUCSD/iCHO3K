function modelFormatted = format_mouse_model(model)

% Replace dashes and other special characters with underscores to match
% formatting convention of Human Recon 2 (humanmetabolism.org)

modelFormatted = model;
mets = model.mets;

mets = regexprep(mets, '-', '_');

modelFormatted.mets = mets;