function [modelData] = getModelData(expressionData,model)

% USAGE: 
% % [modelData] = getModelData(expressionData,model)

% INPUTS:
% % expressionData: a structure containing the entire expression data
                    % % gene: list of gene names
                    % % valuebyTissue: a matrix whose each column
                        % % contains the expression of genes across a tissue/context and 
                        % % each row contains expression of a gene across
                        % % all tissues/contexts.
                    % % genesymbols (optional): any alternative gene names that the
                        % % user may want to keep track of.
                    % % Tissue: names of conditions
% % model:          a genome scale model to be used for integration

% OUTPUTS:
% % modelData:      a structure containing information about the expression of genes in model

% AUTHOR:
% % Chintan Joshi:  for StanDep paper (May 2018)

% if gene ids are in double, convert them to string
if isnumeric(expressionData.gene)
    expressionData.gene = cellstr(char(num2str(expressionData.gene)));
    expressionData.gene = strrep(expressionData.gene,' ','');
end
geneNames = expressionData.gene; 

% construct model data structure
if isfield(expressionData,'symbol')
    modelIdx = ismember(expressionData.gene,model.genes) | ismember(expressionData.symbol,model.genes);
else
    modelIdx = ismember(expressionData.gene,model.genes);
end
modelData.gene = geneNames(modelIdx);
if isfield(expressionData,'meanValue')
    modelData.meanValue = expressionData.meanValue(modelIdx);
end
if isfield(expressionData,'unit')
    modelData.unit = expressionData.unit;
end
if isfield(expressionData,'Tissue')
    modelData.Tissue = expressionData.Tissue;
end
if isfield(expressionData,'genesymbol')
    modelData.genesymbol = expressionData.genesymbol(modelIdx);
end
for j=1:size(expressionData.valuebyTissue,2)
    modelData.value(:,j) = expressionData.valuebyTissue(modelIdx,j);
end

% sore missing gene IDs
modelData.ID_geneMissing = model.genes(~ismember(model.genes,expressionData.gene));
modelData.ID_genePresent = model.genes(ismember(model.genes,expressionData.gene));

% resolve non-uniques
geneids = unique(modelData.gene);
data = zeros(length(geneids),length(modelData.Tissue));
for i=1:length(geneids)
    ix = ismember(modelData.gene,geneids{i});
    if sum(ix)>1
        data(i,:) = mean(modelData.value(ix,:),1);
    else
        data(i,:) = modelData.value(ix,:);
    end
end
modelData.value = data;
modelData.gene = geneids;