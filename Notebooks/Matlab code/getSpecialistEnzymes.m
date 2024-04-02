function [spec] = getSpecialistEnzymes(model)

% USAGE:
% % [spec] = getSpecialistEnzymes(model)

% INPUTS:
% % model:      a COBRA model for which a list of specialist enzymes are desired

% OUTPUTS:
% % spec:       a structure containing information about specialist enzymes

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

% parse and arrange GPRs (needed COBRA toolbox)
parsedGPR = GPRparser(model);
[parsedGPR,ix] = linearization_index(parsedGPR,'rows');
corrRxns = model.rxns(ix);
corrSys = model.subSystems(ix);
ix = find(cellfun(@isempty,parsedGPR));
parsedGPR(ix) = []; corrRxns(ix) = [];
spec.parsedGPR = parsedGPR;
spec.corrRxns = corrRxns;

for i=1:length(parsedGPR)
    if length(parsedGPR)~=1
        parsedGPR{i,1} = strjoin(parsedGPR{i},' & ');
    end
end

% find unique parsedGPR
ugprs = unique(parsedGPR);

% find all other rxns this GPR participates in
rxns = cell(length(ugprs),1);
for i=1:length(ugprs)
    rxns{i,1} = corrRxns(ismember(parsedGPR,ugprs{i}));
    subSystems{i,1} = corrSys(ismember(parsedGPR,ugprs{i}));
end

spec.enzymes = ugprs;
spec.rxns = rxns;
spec.subSystems = subSystems;
nrxns = cellfun(@length,rxns);
spec.enzymes(nrxns>1) = [];
spec.rxns(nrxns>1) = [];
spec.subSystems(nrxns>1) = [];
spec.rxns = linearization_index(spec.rxns,'cols');