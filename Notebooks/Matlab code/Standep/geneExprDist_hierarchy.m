function [clustObj,Z,Hg,Tg,outpermg] = geneExprDist_hierarchy(exprData,removeObjects,edgeX,k,distMethod,linkageMethod)
% USAGE: 
% % [clustObj,Z,Hg,Tg,outpermg,tisMat,rowNames] = geneExprDist_hierarchy(exprData,removeObjects,edgeX,k,distMethod,linkageMethod)

% INPUTS:
% % exprData:       a structure containing expression data (gene/enzyme)
                    % % gene/enzyme: list of gene/enzyme names
                    % % value/valuebyTissue: a matrix whose each column
                        % % contains the expression of genes/enzymes across a tissue/context and 
                        % % each row contains expression of a gene across
                        % % all tissues/contexts.
                    % % rxns:  only needed if enzyme information was used, list of reactions 
                        % % catalyzed by each of the enzymes.
                    % % genesymbols: any alternative gene names that the
                        % % user may want to keep track of.
% % removeObjects:  names of genes or enzymes that need to be removed or
                    % % else give [].
% % edgeX:          give bins acros which distributions are to be
                    % % calculated (required)
% % k:              number of clusters for hierarchical clustering
% % distMethod:     distance method see MATLAB pdist function
% % linkageMethod:  linkage method, see MATLAB linkage function

% OUTPUTS:
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
% % Z:              output of linkage.
% % Hg/Tg/outpermg: outputs of dendrogram (see MATLAB dendrogram function).

% NOTE:
% % The code uses cubehelix for getting colors, please feel free to
% % download it from MATLAB FileExchange (link below):
% % https://www.mathworks.com/matlabcentral/fileexchange/43700-cubehelix-colormap-generator-beautiful-and-versatile

% AUTHOR:
% % Chintan Joshi:  for StanDep paper (May 2018)


if isfield(exprData,'value')
    exprData.valuebyTissue = exprData.value;
    clear exprData.value
end

if ~isempty(removeObjects)
    remIndex = find(ismember(exprData.gene,removeObjects));
    exprData.gene(remIndex) = [];
    exprData.valuebyTissue(remIndex,:) = [];
end
if  isfield(exprData,'gene')
    exprData.gene(sum(exprData.valuebyTissue,2)==0) = [];
elseif isfield(exprData,'enzyme')
    enz=1;
    exprData.enzyme(sum(exprData.valuebyTissue,2)==0) = [];
    exprData.rxns(sum(exprData.valuebyTissue,2)==0) = [];
end
exprData.valuebyTissue(sum(exprData.valuebyTissue,2)==0,:) = [];
nGenes = size(exprData.valuebyTissue,1);
nTis = size(exprData.valuebyTissue,2);

v = log10(exprData.valuebyTissue);

% figure out number of number of columns and rows in figure
nrow = 2;
ncol = 2;
while nrow*ncol <= k
    ncol = ncol + 1;
    if nrow*ncol <= k
        nrow = nrow + 1;
    end
end

Yhist = zeros(nGenes,length(edgeX)-1);
for i=1:nGenes
    temp = v(i,:);
    temp(temp==-inf) = [];
    [Yhist(i,:)] = histcounts(temp,edgeX);
end
Yhist = Yhist/nTis; % normalize by total number of tissues

% hierarchical clustering
Y = pdist(Yhist,distMethod);
Z = linkage(Y,linkageMethod);
leafOrderg = optimalleaforder(Z,Y);
fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
    linkageMethod,distMethod,cophenet(Z,Y));
figure;
[Hg,Tg,outpermg] = dendrogram(Z,k,'Orientation','left','Reorder',leafOrderg);
idxg = Tg;

Cg = zeros(k,length(edgeX)-1);
for i=1:k
    Cg(i,:) = mean(Yhist(idxg==i,:),1);
end
clustObj.Distribution = Yhist;
numObjs = roundn(sum(Yhist,2).*nTis,0);
clustObj.numObjs = numObjs;
if  isfield(exprData,'gene')
    clustObj.objects = exprData.gene;
elseif isfield(exprData,'enzyme')
    clustObj.objects = exprData.enzyme;
    clustObj.objectMaps = exprData.rxns;
end
if isfield(exprData,'genesymbol')
    clustObj.altObjects = exprData.genesymbol;
end
clustObj.C = Cg;
clustObj.cindex = idxg;
clustObj.Data = v;
[idxg,ig] = sort(idxg,'ascend');
Yhist = Yhist(ig,:); 

xtiklab = cell(length(edgeX)-1,1);
for i=1:length(edgeX)-1
    xtiklab{i} = ['[' num2str(edgeX(i)) ',' num2str(edgeX(i+1)) ']'];
end

figure;
cbg = cubehelix(nTis+1,0.6,-0.9,1.9,1.1,[0.2 0.8],[0.4 0.8]);
for i=1:k
    subplot(nrow,ncol,i);
    gid = idxg==i;
    h = plot(1:1:length(xtiklab),Yhist(gid,:)'.*nTis); hold on
    ax = gca; ax.TickDir = 'out';
    set(h,{'Color'},num2cell(cbg(nTis+1-numObjs(gid),:),2));
    if i==k
        colormap(cubehelix([],0.6,-0.9,1.9,1.1,[0.2 0.8],[0.4 0.8]));
        cbi = colorbar('Position',[0.925 0.18 0.02 0.7]);
        cbi.Ticks = linspace(0,1,5);
        cbi.TickLabels = num2cell(nTis:-nTis/4:0);
        cbi.Label.String = 'Number of tissues gene is expressed in';
        cbi.Label.FontSize = 14;
        cbi.TickDirection = 'out';
    end
    plot(1:1:length(xtiklab),Cg(i,:)*nTis,'Color','k','LineWidth',2); hold off
    xlim([1 (length(edgeX)-1)]);
    xticks(1:1:(length(edgeX)-1));
    xticklabels(xtiklab);
    xlabel('Bins of log_{10}[expr.]');
    xtickangle(45);
    clustObj.numObjInClust(i,1) = sum(gid);
    if isfield(exprData,'gene')
        title(['Cluster # ' num2str(i) ' (' num2str(sum(gid)) ')']);
    elseif isfield(exprData,'enzyme')
        title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' enzymes']);
    end
    ylim([0 nTis]);
    yticks(0:0.25*nTis:nTis);
    yticklabels(0:0.25*nTis:nTis);
    ylabel('# of tissues');
end

figure;
tisMat = []; % initialize tissueMatrix
rowNames = [];
r = mod(nTis,4);
for i=1:k
    gid = clustObj.cindex==i;
    ngenes = sum(gid);
    geneTisMat = zeros(nTis,length(xtiklab));
    for j=1:nTis
        temp = v(gid,:);
        temp = temp(:,j);
        temp(temp==-inf) = [];
        geneTisMat(j,:) = histcounts(temp,edgeX)';
    end;
    tisMat = [tisMat;geneTisMat/ngenes]; % join tissueMatrix
    myStrings = strcat(exprData.Tissue,'_',num2str(i));
    rowNames = [rowNames;myStrings];
    geneTisMat = 1-geneTisMat/ngenes;
    subplot(nrow,ncol,i);
    imagesc(geneTisMat);
    caxis([0 1]);
    ax = gca; ax.TickDir = 'out';
    ylim([0 nTis]);
    yticks(0:0.25*(nTis-r):(nTis-r));
    yticklabels(0:0.25*(nTis-r):(nTis-r));
    ylabel('cell type #');
    xlim([0.5 (length(edgeX)-0.5)]);
    xticks(1:1:(length(edgeX)-1));
    xticklabels(xtiklab);
    xlabel('Bins of log_{10}[expr.]');
    xtickangle(45);
    colormap(cubehelix([],0.1,0.4,2,1,[0.1 1],[0 1]));
    if isfield(exprData,'gene')
        title(['Cluster # ' num2str(i) ' (' num2str(sum(gid)) ')']);
    elseif isfield(exprData,'enzyme')
        title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' enzymes']);
    end
end
colormap(cubehelix([],0.1,0.4,2,1,[0.1 1],[0 1]));
cbi = colorbar('Position',[0.925 0.18 0.02 0.7]);
cbi.Ticks = linspace(0,1,6);
cbi.TickLabels = num2cell(1:-1/5:0);
cbi.Label.String = 'fraction of genes in the cluster';
cbi.Label.FontSize = 14;
cbi.TickDirection = 'out';