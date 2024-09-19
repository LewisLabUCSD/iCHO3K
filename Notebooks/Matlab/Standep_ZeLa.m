warning('off', 'all');

% Initialize COBRA
initCobraToolbox(0)

path = '/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells';
addpath('/Users/pablodigiusto/Documents/GitHub/Whole-Cell-Network-Reconstruction-for-CHO-cells/Notebooks/Matlab/Standep');

% Load iCHO3644 model
model = readCbModel(fullfile(path, 'Notebooks', 'iCHO3595_unblocked.mat'));

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
    'S1', 'WT_P2_Bio141'; 'S7', 'WT_P2_Bio142'; 'S13', 'WT_P2_Bio143';
    'S2', 'WT_P4_Bio141'; 'S8', 'WT_P4_Bio142'; 'S14', 'WT_P4_Bio143';
    'S3', 'WT_P6_Bio141'; 'S9', 'WT_P6_Bio142'; 'S15', 'WT_P6_Bio143';
    'S4', 'WT_P8_Bio141'; 'S10', 'WT_P8_Bio142'; 'S16', 'WT_P8_Bio143';
    'S5', 'WT_P12_Bio141';'S11', 'WT_P12_Bio142';
    'S6', 'WT_P14_Bio141';'S12', 'WT_P14_Bio142';'S17', 'WT_P14_Bio143';
    'S18', 'ZeLa_P4_Bio144'; 'S23', 'ZeLa_P4_Bio145'; 'S28', 'ZeLa_P4_Bio146'; 'S34', 'ZeLa_P4_Bio147'; 'S39', 'ZeLa_P4_Bio148';
    'S19', 'ZeLa_P6_Bio144'; 'S29', 'ZeLa_P6_Bio146'; 'S35', 'ZeLa_P6_Bio147'; 'S40', 'ZeLa_P6_Bio148';
    'S20', 'ZeLa_P8_Bio144'; 'S24', 'ZeLa_P8_Bio145'; 'S30', 'ZeLa_P8_Bio146'; 'S36', 'ZeLa_P8_Bio147'; 'S41', 'ZeLa_P8_Bio148';
    'S25', 'ZeLa_P12_Bio145';'S31', 'ZeLa_P12_Bio146';'S42', 'ZeLa_P12_Bio148';
    'S21', 'ZeLa_P14_Bio144';'S26', 'ZeLa_P14_Bio145';'S32', 'ZeLa_P14_Bio146';'S37', 'ZeLa_P14_Bio147';'S43', 'ZeLa_P14_Bio148';
    'S22', 'ZeLa_P2_Bio145'; 'S27', 'ZeLa_P2_Bio146'; 'S33', 'ZeLa_P2_Bio147'; 'S38', 'ZeLa_P2_Bio148';
};

% Prepare the expression data matrix and cell names
expressionDataMatrix = table2array(dataTable(:, 2:end)); % Assuming the first column is gene names
cellNames = dataTable.Properties.VariableNames(2:end);

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


% extract expression data of the genes in the model
modelData = getModelData(expressionData,model);


% Calculate enzymes in the model
spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);

% Calculate enzyme expression
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

edgeX = [-2 -1 0 1 2 2.5 3 4]; % bins  
k = 3; % what value should we use
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering

% Calculate clusters of enzyme expression
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);

% calculate active reaction lists as binary matrix
coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],true,0,[1 1]);

% calculate Jaccard similarity between core reaction lists
calcJaccardSimilarity(coreRxnMat,enzymeData.Tissue,'matrix',true);

% Calculate ubiquity scores for mCADRE
[ubiScores, uScore] = getUbiquityScore(clustObj, edgeX, model); % calculate ubiquity score

% Create the ubiData structure
ubiData = struct;
ubiData.ubiScores = ubiScores;
ubiData.uScore = uScore;
ubiData.rxns = model.rxns;
ubiData.Condition = cellNames;

% Save the outputs to a .mat file
output_path = fullfile(path, 'Data/Context_specific_models', 'ubiData.mat');
save(output_path, 'ubiData');

