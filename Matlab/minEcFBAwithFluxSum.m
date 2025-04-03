function [modelIrrev, Solutions,fluxSum] = minEcFBAwithFluxSum(model, biomass, sample_size)
% INPUT
%  model             COBRA model structure with 3 additional vectors of same size as 'rxns':
%                    kcat_f, kcat_b, molwt (if any of the value unknown,
%                    provide '0')
%                    kcat units should be '1/s' and mol wt in 'Dalton'
%  frac              fraction of enzymatic mass in overall dry cell weight   
%  biomass           name of biomass reaction (to be excluded from enzyme
%                    capacity flux constraint)
% 
% OUTPUT
%  modelIrrev    Model in irreversible format with enzyme constraint added
%                as pseudo reaction
%  solution
%    f         Objective value
%    x         Primal
%    y         Dual
%    w         Reduced costs
%    s         Slacks
%    stat      Solver status in standardized form
%               1   Optimal solution
%               2   Unbounded solution
%               0   Infeasible
%              -1  No solution reported (timelimit, numerical problem etc)
% Lokanand Koduru            10/03/18
% Meiyappan Lakshmanan       10/04/18 Generalized the code

%% convert to irreversible format
[modelIrrev,~,~,irrev2rev] = convertToIrreversible_local(model);
modelIrrev.description='model_with_EnzCon';
revFlag='false';
%% identify reactions with kcat values
kcat_rxns = model.rxns(sum([model.kcat_f,model.kcat_b], 2)~=0);
%% identify reactions in Irreversible model with kcat values
kcat_rxns_irrev = modelIrrev.rxns(ismember(irrev2rev,find(ismember(model.rxns,kcat_rxns))));
%% assign kcat values to reactions in Irreversible model with kcat values
kcat_irrev = zeros(length(kcat_rxns_irrev),1);
for i=1:1:length(kcat_rxns_irrev)
    if model.rev(find(ismember(model.rxns,kcat_rxns_irrev(i)))) == 0
        kcat_irrev(i) = model.kcat_f(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i)))));
    elseif model.kcat_b(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i))))) ~= 0
        kcat_irrev(i) = model.kcat_b(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i)))));
    else
        kcat_irrev(i) = model.kcat_f(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i)))));
    end
end
%% assign mw values to reactions in Irreversible model with kcat values
for i=1:1:length(kcat_rxns_irrev)
    mw_irrev(i) = model.molwt(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i)))));
end
kcat_mw = zeros(length(kcat_rxns_irrev),1);
%% calculate kcat/mw values
for i=1:1:length(kcat_rxns_irrev)
    kcat_mw(i) = (mw_irrev(i)/1000)/(kcat_irrev(i)*3600);
end

kcat_rxns_irrev = kcat_rxns_irrev(~isinf(kcat_mw));
kcat_mw = kcat_mw(~isinf(kcat_mw));
kcat_rxns_irrev = kcat_rxns_irrev(kcat_mw~=0);
kcat_mw = kcat_mw(kcat_mw~=0);
kcat_rxns_irrev = kcat_rxns_irrev(~isnan(kcat_mw));
kcat_mw = kcat_mw(~isinf(kcat_mw));
%% add Dummy reaction to model repreenting the enzyme constraint
lowerBound=0;
upperBound=1000;
modelIrrev=addReaction(modelIrrev,'EnzCon_C_Rxn',{'EnzCon_Met[c]'},-1,revFlag,lowerBound,upperBound,0,'','','','');
modelIrrev = changeObjective(modelIrrev,'EnzCon_C_Rxn');
EnzCon_MetInd=find(ismember(modelIrrev.mets,'EnzCon_Met[c]'));
Rxn_WithEnzConInd=find(ismember(modelIrrev.rxns,kcat_rxns_irrev));
selExc=findExcRxns(modelIrrev);
ExRxnInd=find(selExc);
BiomassRxnInd=find(ismember(modelIrrev.rxns,biomass));
Solutions=zeros(length(modelIrrev.rxns),sample_size);
fluxSum = zeros(length(modelIrrev.mets),sample_size);           %Added by VIKASH & KANNAN
EnzConSamples=zeros(length(modelIrrev.rxns)-1,sample_size);
RandomEnzConCoeff=zeros(length(modelIrrev.rxns),sample_size);
sort_kcat_mw = sortrows(kcat_mw);
for i=1:1:sample_size
    RandomEnzConCoeff(:,i)=randsample(kcat_mw,length(modelIrrev.rxns),true);
end
for j=1:1:sample_size
%     changeCobraSolver('gurobi');
    modelIrrev1 = modelIrrev;
    kcat_mw1 = kcat_mw;
    modelIrrev1.S(EnzCon_MetInd,:)=transpose(RandomEnzConCoeff(:,j));
    %% Removing EnzCon coeff for biomass equation
    modelIrrev1.S(EnzCon_MetInd,BiomassRxnInd)=0;

    %% Replacing the randomly assigned EnzCon coeff with the original coefficients for those reactions whose EnzCon coeffs are available (calculated list) 
    for k=1:length(Rxn_WithEnzConInd)
            modelIrrev1.S(EnzCon_MetInd,Rxn_WithEnzConInd(k))=kcat_mw1(find(ismember(kcat_rxns_irrev,modelIrrev1.rxns(Rxn_WithEnzConInd(k)))));
    end
    %% Removing EnzCon coeff for exchange reactions
    modelIrrev1.S(EnzCon_MetInd,ExRxnInd)=0;
    modelIrrev1.S(EnzCon_MetInd,find(ismember(modelIrrev1.rxns,'EnzCon_C_Rxn')))=-1;
    solution=optimizeCbModel(modelIrrev1,'min');
    if(isempty(solution.x))==0
        Solutions(:,j)=solution.x;
       fs = sum(abs((modelIrrev1.S)*(diag(solution.x)))' * 0.5)';           %Added by VIKASH & KANNAN
        fluxSum(:,j)=fs;
    end
end
valid_sol = find(Solutions(BiomassRxnInd,:)>0);
Solutions = Solutions(:,valid_sol);
end

function [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible_local(model)
%convertToIrreversible Convert model to irreversible format
%
% [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model)
%
%INPUT
% model         COBRA model structure
%
%OUTPUTS
% modelIrrev    Model in irreversible format
% matchRev      Matching of forward and backward reactions of a reversible
%               reaction
% rev2irrev     Matching from reversible to irreversible reactions
% irrev2rev     Matching from irreversible to reversible reactions
%
% Uses the reversible list to construct a new model with reversible
% reactions separated into forward and backward reactions.  Separated
% reactions are appended with '_f' and '_b' and the reversible list tracks
% these changes with a '1' corresponding to separated forward reactions.
% Reactions entirely in the negative direction will be reversed and
% appended with '_r'.
%
% written by Gregory Hannum 7/9/05
%
% Modified by Markus Herrgard 7/25/05
% Modified by Jan Schellenberger 9/9/09 for speed.

%declare variables
modelIrrev.S = spalloc(size(model.S,1),0,2*nnz(model.S));
modelIrrev.rxns = [];
modelIrrev.rev = zeros(2*length(model.rxns),1);
modelIrrev.lb = zeros(2*length(model.rxns),1);
modelIrrev.ub = zeros(2*length(model.rxns),1);
modelIrrev.c = zeros(2*length(model.rxns),1);
matchRev = zeros(2*length(model.rxns),1);

nRxns = size(model.S,2);
irrev2rev = zeros(2*length(model.rxns),1);

%loop through each column/rxn in the S matrix building the irreversible
%model
cnt = 0;
for i = 1:nRxns
    cnt = cnt + 1;
    
    %expand the new model (same for both irrev & rev rxns  
    modelIrrev.rev(cnt) = model.rev(i);
    irrev2rev(cnt) = i;

    % Check if reaction is declared as irreversible, but bounds suggest
    % reversible (i.e., having both positive and negative bounds
    if (model.ub(i) > 0 && model.lb(i) < 0) && model.rev(i) == false
        model.rev(i) = true;
        warning(cat(2,'Reaction: ',model.rxns{i},' is classified as irreversible, but bounds are positive and negative!'))

    end
   
    % Reaction entirely in the negative direction
    if (model.ub(i) <= 0 && model.lb(i) < 0)
        % Retain original bounds but reversed
        modelIrrev.ub(cnt) = -model.lb(i);
        modelIrrev.lb(cnt) = -model.ub(i);
        % Reverse sign
        modelIrrev.S(:,cnt) = -model.S(:,i);
        modelIrrev.c(cnt) = -model.c(i);
        modelIrrev.rxns{cnt} = [model.rxns{i} '_r'];
        model.rev(i) = false;
        modelIrrev.rev(cnt) = false;
    else
        % Keep positive upper bound
        modelIrrev.ub(cnt) = model.ub(i);
        %if the lb is less than zero, set the forward rxn lb to zero 
        if model.lb(i) < 0
            modelIrrev.lb(cnt) = 0;
        else
            modelIrrev.lb(cnt) = model.lb(i);
        end
        modelIrrev.S(:,cnt) = model.S(:,i);
        modelIrrev.c(cnt) = model.c(i);
        modelIrrev.rxns{cnt} = model.rxns{i};

    end

   
    %if the reaction is reversible, add a new rxn to the irrev model and
    %update the names of the reactions with '_f' and '_b'
    if model.rev(i) == true
        cnt = cnt + 1;
        matchRev(cnt) = cnt - 1;
        matchRev(cnt-1) = cnt;
        modelIrrev.rxns{cnt-1} = [model.rxns{i} '_f'];
        modelIrrev.S(:,cnt) = -model.S(:,i);
        modelIrrev.rxns{cnt} = [model.rxns{i} '_b'];
        modelIrrev.rev(cnt) = true;
        if model.lb(i) >0 && model.ub(i)>0
            modelIrrev.lb(cnt-1) = model.lb(i);
            modelIrrev.ub(cnt-1) = model.ub(i);
            modelIrrev.lb(cnt) = 0;
            modelIrrev.ub(cnt) = 0;
        else
            modelIrrev.lb(cnt) = 0;
            modelIrrev.ub(cnt) = -model.lb(i);
        end
        modelIrrev.c(cnt) = 0;
        rev2irrev{i} = [cnt-1 cnt];
        irrev2rev(cnt) = i;
    else
        matchRev(cnt) = 0;
        rev2irrev{i} = cnt;
    end
end

rev2irrev = columnVector(rev2irrev);
irrev2rev = irrev2rev(1:cnt);
irrev2rev = columnVector(irrev2rev);

% Build final structure
modelIrrev.S = modelIrrev.S(:,1:cnt);
modelIrrev.ub = columnVector(modelIrrev.ub(1:cnt));
modelIrrev.lb = columnVector(modelIrrev.lb(1:cnt));
modelIrrev.c = columnVector(modelIrrev.c(1:cnt));
modelIrrev.rev = modelIrrev.rev(1:cnt);
modelIrrev.rev = columnVector(modelIrrev.rev == 1);
modelIrrev.rxns = columnVector(modelIrrev.rxns); 
modelIrrev.mets = model.mets;
matchRev = columnVector(matchRev(1:cnt));
modelIrrev.match = matchRev;
if (isfield(model,'b'))
    modelIrrev.b = model.b;
end
if isfield(model,'description')
    modelIrrev.description = [model.description ' irreversible'];
end
if isfield(model,'subSystems')
    modelIrrev.subSystems = model.subSystems(irrev2rev);
end
if isfield(model,'genes')
    modelIrrev.genes = model.genes;
    genemtxtranspose = model.rxnGeneMat';
    modelIrrev.rxnGeneMat = genemtxtranspose(:,irrev2rev)';
    modelIrrev.rules = model.rules(irrev2rev);
end
modelIrrev.reversibleModel = false;

end