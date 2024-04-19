function [GPRrxns, GPRmat] = parse_gprs(model)

numGPRs = cellfun('length', regexp(model.grRules, 'or'));

noRules = cellfun('isempty', model.grRules);
numGPRs(~noRules) = numGPRs(~noRules) + 1;

numRows = sum(numGPRs); numCols = max(numGPRs);

GPRmat = repmat({''}, numRows, numCols);
GPRrxns = repmat({''}, numRows, 1);

count = 1;
for i = 1:numel(model.grRules)
    if numel(model.grRules{i})
        ruleGPRs = regexp(model.grRules{i}, 'or', 'split');
        ruleGPRs = regexprep(ruleGPRs, '[\s\(\)]', '');
        for j = 1:numel(ruleGPRs)
            GPR = regexp(ruleGPRs{j}, 'and', 'split');
            GPRmat(count,1:numel(GPR)) = GPR;
            GPRrxns(count) = model.rxns(i);
            count = count + 1;
        end
    end
end

GPRmat(:,sum(~cellfun('isempty', GPRmat), 1) == 0) = [];
GPRmat = regexprep(GPRmat, '\.[0-9]', '');
