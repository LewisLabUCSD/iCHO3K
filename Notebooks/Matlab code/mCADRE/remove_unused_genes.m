function modelNew = remove_unused_genes(model)

%Updating genes after reaction removal process finished
genes = unique(model.genes);

modelNew = model;

modelNew.rules = [];
modelNew.rxnGeneMat = [];
modelNew.genes = [];

for m = 1:numel(modelNew.rxns)
    [~, modelNew] = evalc(['changeGeneAssociation(modelNew, ', ...
        'modelNew.rxns{m, 1}, modelNew.grRules{m, 1}, genes, genes);']);
end
