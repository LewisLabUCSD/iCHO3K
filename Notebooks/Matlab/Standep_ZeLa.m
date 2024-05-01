warning('off', 'all');

% Initialize COBRA
initCobraToolbox

path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';
addpath('/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells/Notebooks/Matlab/Standep');

% Load iCHO3644 model
model = readCbModel(fullfile(path, 'Notebooks', 'iCHO3644_unblocked.mat'));

% Create a dictionary with the gene names of our model as keys and the gene IDs as values

genes = model.genes; % or model.rxns or another field depending on your model
geneNames = model.geneNames; % Adjust field name as per your model's structure

% Initialize a map (dictionary) to hold gene names as keys and genes IDs as values
geneDict = containers.Map('KeyType', 'char', 'ValueType', 'char');

% Iterate over the genes
for i = 1:length(genes)
    geneID = genes{i};
    geneName = geneNames{i};
    geneDict(geneName) = geneID; % Use gene name as key, gene ID as value
end

% Loading the RNA-seq data
rnaseq_path = fullfile(path, 'Data/Zela Data', '20200307_Bio141-148_merged.tpm.tsv');
dataTable = readtable(rnaseq_path, 'FileType', 'text', 'Delimiter', '\t');

% Extracting gene names if the first column
geneNames = dataTable{:, 1};

% Map gene IDs from the model to the RNA-Seq data
for i = 1:length(geneNames)
    currentGeneName = geneNames{i};
    if isKey(geneDict, currentGeneName)
        % If the current gene name exists in the dictionary, replace it with the gene ID
        updatedGeneNames{i} = geneDict(currentGeneName);
    else
        % If the gene name does not exist in the dictionary, keep the original name
        updatedGeneNames{i} = currentGeneName;
    end
end

% Sample Identifiers and Their Conditions
sampleConditions = {
    'S1', 'WT_P2'; 'S7', 'WT_P2'; 'S13', 'WT_P2';
    'S2', 'WT_P4'; 'S8', 'WT_P4'; 'S14', 'WT_P4';
    'S3', 'WT_P6'; 'S9', 'WT_P6'; 'S15', 'WT_P6';
    'S4', 'WT_P8'; 'S10', 'WT_P8'; 'S16', 'WT_P8';
    'S5', 'WT_P12';'S11', 'WT_P12';
    'S6', 'WT_P14';'S12', 'WT_P14';'S17', 'WT_P14';
    'S18', 'ZeLa_P4'; 'S23', 'ZeLa_P4'; 'S28', 'ZeLa_P4'; 'S34', 'ZeLa_P4'; 'S39', 'ZeLa_P4';
    'S19', 'ZeLa_P6'; 'S29', 'ZeLa_P6'; 'S35', 'ZeLa_P6'; 'S40', 'ZeLa_P6';
    'S20', 'ZeLa_P8'; 'S24', 'ZeLa_P8'; 'S30', 'ZeLa_P8'; 'S36', 'ZeLa_P8'; 'S41', 'ZeLa_P8';
    'S25', 'ZeLa_P12';'S31', 'ZeLa_P12';'S42', 'ZeLa_P12';
    'S21', 'ZeLa_P14';'S26', 'ZeLa_P14';'S32', 'ZeLa_P14';'S37', 'ZeLa_P14';'S43', 'ZeLa_P14';
    'S22', 'ZeLa_P2'; 'S27', 'ZeLa_P2'; 'S33', 'ZeLa_P2'; 'S38', 'ZeLa_P2';
};

% Prepare the expression data matrix and cell names
expressionDataMatrix = table2array(dataTable(:, 2:end)); % Assuming the first column is gene names
cellNames = dataTable.Properties.VariableNames(2:end); % Adjust if the structure is different

% Adjust cellNames based on sampleConditions
for i = 1:size(sampleConditions, 1)
    idx = find(strcmp(cellNames, sampleConditions{i, 1}));
    if ~isempty(idx)
        cellNames{idx} = sampleConditions{i, 2};
    end
end

% Creating the expressionData Structure
expressionData = struct;
expressionData.gene = updatedGeneNames;
expressionData.valuebyTissue = expressionDataMatrix;
expressionData.Tissue = cellNames;


% Get the mean values of gene expression per condition

% Find unique condition names
uniqueConditions = unique(cellNames);

% Initialize a matrix to store the mean values for each condition
meanExpressionDataMatrix = zeros(size(expressionDataMatrix, 1), length(uniqueConditions));

% Loop over each unique condition to calculate the mean expression
for i = 1:length(uniqueConditions)
    % Find columns that belong to the current condition
    conditionCols = strcmp(cellNames, uniqueConditions{i});
    
    % Calculate the mean across these columns for each gene
    meanExpressionDataMatrix(:, i) = mean(expressionDataMatrix(:, conditionCols), 2);
end

% Creating the meanexpressionData Structure
meanexpressionData = struct;
meanexpressionData.gene = updatedGeneNames;
meanexpressionData.valuebyTissue = meanExpressionDataMatrix;
meanexpressionData.Tissue = uniqueConditions;


% extract expression data of the genes in the model
modelData = getModelData(meanexpressionData,model);


% Calculate enzymes in the model
spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);

% Calculate enzyme expression
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

edgeX = [-2 -1 0 1 2 2.5 3 4]; % bins  
k = 20; % what value should we use
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering

% Calculate clusters of enzyme expression
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);

% calculate active reaction lists as binary matrix
coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],true,0,[1 1]);

% calculate Jaccard similarity between core reaction lists
calcJaccardSimilarity(coreRxnMat,enzymeData.Tissue,'matrix',true);

% Calculate ubiquity scores for mCADRE
[ubiScore,uScore] = getUbiquityScore(clustObj,edgeX,model); % calculate ubiquity score

size(ubiScore) % ubiquity score for reactions

size(uScore) % this is in theory ubiquity scores for enzymes.. ??

% Save the outputs to a .mat file
output_path = fullfile(path, 'Data/Context_specific_models/UbiquityScores.mat');
save(output_path, 'ubiScore', 'uScore');

