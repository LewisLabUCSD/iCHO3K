function modelFormatted = format_human_model(model)

% Replace dashes and other special characters with underscores to match
% formatting convention of Human Recon 2 (humanmetabolism.org)

modelFormatted = model;
mets = model.mets;

mets = regexprep(mets, '-', '_');
mets = regexprep(mets, '_DASH_', '_'); 
mets = regexprep(mets, '_FSLASH_', '_'); 
mets = regexprep(mets,' _COMMA_', '_');

modelFormatted.mets = mets;
