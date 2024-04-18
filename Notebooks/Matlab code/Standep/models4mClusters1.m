function [geneTis,enzTis,cutOff,thr,enzSel,tisClust] = models4mClusters1(clustObj,tisNames,model,xedge,folderName,cutOff,figFlag,adjustValue,weightage)
% USAGE:
% % [geneTis,enzTis,cutOff,thr,enzSel,tisClust] = models4mClusters1(clustObj,tisNames,model,xedge,folderName,cutOff,figFlag,adjustValue,weightage)

% INPUTS:
% % clustObj:       a structure containing information about the clusters
                    % % Distribution: contains gene-specific distributions using bins provided
                    % % objects: gene/enzyme names provided
                    % % objectMaps: only exists if enzymes were used, contains list of
                        % % reactions catalyzed by the model. This is also must be calculated before
                        % % and provided in exprData.
                    % % cindex: cluster indices [1,k]
                    % % Data: Original data
                    % % C: mean vector of the cluster
                    % % numObjInClust: number of genes/enzymes in each cluster
                    % % altObjects: any other alternative gene names also supplied with the exprData
% % tisNames:       list of tissue or context names/ids
% % model:          genome-scale model to be used
% % xedge:          the bins of expression value to be used generating figure
% % foldername:     if not writing a file, []; else, give folder location (even if using current working directory)
% % cutOff:         if using StanDep code for calculating thresholds, [].
% % figFlag:        print the figure, true; else, false
% % adjustValue:    if want to adjust thresholds across all clusters by a
                        % % certain value provide it here, else 0.
% % weightage:      to give different weightage to standard deviation and
                        % % mean terms, else say [1 1].
                        
% OUTPUTS:
% % geneTis:        if using enzyme clusters in clustObj, contains
                        % % selection of reactions; if using genes in clustObj, contains selection of
                        % % genes across tissues.
% % enzTis:         contains enzymes used in each tissue, ignore this if
                        % % using gene clusters
% % cutOff:         cutOffs calculated using clusterVariability1.m in percentile of each cluster
% % thr:            threshold in absolute values of expression units
% % enzSel:         the enzyme in the tissue comes from which cluster
% % tisClust:       a tissue cluster object structure
                    % % outperm/H/T:outputs of dendrogram (see MATLAB function)
                    % % tisNames:   tissue names

% AUTHOR:
% % Chintan Joshi:  for StanDep paper (May 2018)

muData = clustObj.Data; muData(muData==-inf) = [];
muData = reshape(muData,numel(muData),1);
muData = mean(muData);
if isempty(cutOff)
    [cutOff,thr] = clusterVariability1(clustObj,xedge,figFlag,adjustValue,weightage);
end
nTis = size(clustObj.Data,2);
nClust = size(clustObj.C,1);
tisNames = strrep(tisNames,'_',' ');
geneTis = false(length(clustObj.objects),nTis);
enzSel = zeros(length(clustObj.objects),nTis);
for i=1:nClust
    gid = clustObj.cindex==i;
    clustMat = false(sum(gid),nTis);
    clustData = clustObj.Data(gid,:);
    clustData = reshape(clustData,numel(clustData),1);
    colgene = repmat(clustObj.objects(gid),nTis,1);
    coltis = repmat(1:1:nTis,sum(gid),1);
    coltis = reshape(coltis,numel(coltis),1);
    remi = clustData==-inf;
    coltis(remi) = []; colgene(remi) = [];
    clustData(remi) = [];
    thr(i,1) = prctile(clustData,roundn(cutOff(i),-2));
    ind = clustData >= thr(i);
    fprintf('fraction selected = %0.4f\n',sum(ind)/(sum(gid)*nTis));
    coltis = coltis(ind);
    colgene = colgene(ind);
    for j=1:nTis
        if sum(coltis==j)~=0
            geneTis(:,j) = geneTis(:,j) | ismember(clustObj.objects,unique(colgene(coltis==j)));
            enzSel(ismember(clustObj.objects,unique(colgene(coltis==j))),j) = i;
            clustMat(:,j) = ismember(clustObj.objects(gid),unique(colgene(coltis==j)));
        end
    end
    if ~isempty(folderName)
        Jc = getJaccardSimMatrix(clustMat);
        writeMyMatrix(Jc,strcat(folderName,'/C',num2str(i),'.txt'));
        if i==1
            Z = linkage(Jc,'complete');
            Y = pdist(Jc,'euclidean');
            leafOrderg = optimalleaforder(Z,Y);
            fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
                'complete','euclidean',cophenet(Z,Y));
            if figFlag
                figure;
                [Hg,Tg,outpermg] = dendrogram(Z,0,'Orientation','left','labels',tisNames,'Reorder',leafOrderg);
                title('Tissues clustered based on metabolic similarity in Cluster 1');
            end
        end
    end
end
enzTis = geneTis;
J = getJaccardSimMatrix(geneTis);
if ~isempty(folderName)
    writeMyMatrix(J,strcat(folderName,'/Call.txt'));
end
tisClust = [];
if figFlag
    Z = linkage(J,'complete');
    Y = pdist(J,'euclidean');
    leafOrderg = optimalleaforder(Z,Y);
    fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
        'complete','euclidean',cophenet(Z,Y));
    figure;
    [Hg,Tg,outpermg] = dendrogram(Z,0,'Orientation','left','labels',tisNames,'Reorder',leafOrderg);
    tisClust.outperm = outpermg; tisClust.H = Hg; tisClust.T = Tg;
    tisClust.names = tisNames;
    title('Tissues clustered based on metabolic similarity');
    outpermg = outpermg(sort(1:1:nTis,'descend'));
    J = J(outpermg,outpermg);
    tisNames = tisNames(outpermg);

    figure;
    imagesc(J);
    ax = gca; ax.TickDir = 'out'; ax.DataAspectRatioMode = 'manual';
    xticks(1:1:nTis); yticks(1:1:nTis);
    xticklabels(tisNames);
    yticklabels(tisNames);
    xtickangle(45);
    if ~isfield(clustObj,'objectMaps')
        title('Jaccard Similarity for genes in tissues');
    else
        title('Jaccard Similarity for enzymes in tissues');
    end
    cbi = colorbar;
    cbi.TickDirection = 'out';
    colormap jet
end

if ~isempty(model)
    if ~isfield(clustObj,'objectMaps')
        model.genes = unique(model.genes);
        activeGenes = cell(1,nTis);
        for i=1:nTis
            activeGenes{i} = find(ismember(model.genes,clustObj.objects(geneTis(:,i))));
        end
        geneTis = activeGenes;
    else
        activeRxns = false(length(model.rxns),nTis);
        for i=1:nTis
            mylist = clustObj.objectMaps(geneTis(:,i));
            l1 = mylist(~cellfun(@iscell,mylist));
            l2 = linearization_index(mylist(cellfun(@iscell,mylist)),'cols');
            mylist = unique([l1; l2]);
            activeRxns(:,i) = ismember(model.rxns,mylist);
        end
        geneTis = activeRxns;
    end
end


function J = getJaccardSimMatrix(A)

n = size(A,2);
J = zeros(n,n);
for i=1:n
    for j=i:n
        J(j,i) = sum(A(:,i) & A(:,j))/sum(A(:,i) | A(:,j));
        J(i,j) = sum(A(:,i) & A(:,j))/sum(A(:,i) | A(:,j));
    end
end

function writeMyMatrix(A,fileName)

fid = fopen(fileName,'wt');

for ii = 1:size(A,1)
    fprintf(fid,'%0.2f\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);