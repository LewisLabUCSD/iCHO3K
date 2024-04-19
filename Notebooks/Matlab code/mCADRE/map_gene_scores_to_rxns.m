function U_GPR = map_gene_scores_to_rxns(model, G, U, GPRmat)

if ~iscellstr(model.genes)
    genes = strtrim(cellstr(num2str(model.genes)));
else genes = model.genes;
end

genes = regexprep(genes, '\.[0-9]', '');
U_GPR = nan(size(GPRmat));
[geneInt, G_idx] = intersect(G, genes);
for i = 1:numel(geneInt)
    gene_GPR = strcmp(GPRmat, geneInt{i});
    U_GPR(gene_GPR) = U(G_idx(i));
end

% Penalize genes with zero expression, such that corresponding reactions
% will be ranked lower than non-gene associated reactions.
U_GPR(U_GPR == 0) = -1e-6;
