function [prom] = getPromEnzymes(model)

% USAGE:
% % [prom] = getPromEnzymes(model)

% INPUTS:
% % model:      a COBRA model for which a list of promiscuous enzymes are desired

% OUTPUTS:
% % prom:       a structure containing information about promiscuous enzymes

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

% parse and arrange GPRs
parsedGPR = GPRparser(model);
[parsedGPR,ix] = linearization_index(parsedGPR,'rows');
corrRxns = model.rxns(ix);
corrSys = model.subSystems(ix);
ix = find(cellfun(@isempty,parsedGPR));
parsedGPR(ix) = []; corrRxns(ix) = [];
prom.parsedGPR = parsedGPR;
prom.corrRxns = corrRxns;

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

prom.enzymes = ugprs;
prom.rxns = rxns;
prom.subSystems = subSystems;
prom.nrxns = cellfun(@length,rxns);
prom.enzymes(prom.nrxns==1) = [];
prom.rxns(prom.nrxns==1) = [];
prom.subSystems(prom.nrxns==1) = [];
prom.nrxns(prom.nrxns==1) = [];