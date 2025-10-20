function [modelIrrev, Solutions,fluxSum] = maxEcFBAwithFluxSum(model, frac, biomass, sample_size)
% function [modelIrrev, Solutions, ShadowPrices, ReducedCosts, EnzConSamples] = enzymeconstrainedfba_cho(model, frac, biomass, sample_size)
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

% altering the bounds of rxnID : 'EnzCon_C_Rxn' [line: 75-76]
%% convert to irreversible format

%modelIrrev = changeObjective(model,biomass);%%edited by kannan 
[modelIrrev,~,~,irrev2rev] = convertToIrreversible_local(model);
% sol = optimizeCbModel(modelIrrev);
% modelIrrev = changeRxnBounds(modelIrrev,'biomass_reaction',0.1*sol.f,'b');
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
%% calculate kcat/mw values %correction: mw/kcat --done by Bharathi
for i=1:1:length(kcat_rxns_irrev)
    kcat_mw(i) = (mw_irrev(i)/1000)/(kcat_irrev(i)*3600);
end

kcat_rxns_irrev = kcat_rxns_irrev(~isinf(kcat_mw));
kcat_mw = kcat_mw(~isinf(kcat_mw));
kcat_mw= kcat_mw(~isnan(kcat_mw)); %added by Trinankur and Bharathi
kcat_rxns_irrev = kcat_rxns_irrev(kcat_mw~=0);
kcat_mw = kcat_mw(kcat_mw~=0);
kcat_rxns_irrev = kcat_rxns_irrev(~isnan(kcat_mw));
kcat_mw = kcat_mw(~isinf(kcat_mw));

%% add Dummy reaction to model repreenting the enzyme constraint
lowerBound=0;
upperBound=frac;
modelIrrev=addReaction(modelIrrev,'EnzCon_C_Rxn',{'EnzCon_Met[c]'},-1,revFlag,lowerBound,upperBound,0,'','','','');
modelIrrev = changeRxnBounds(modelIrrev,'EnzCon_C_Rxn',0.1,'u'); % tuned by Trinankur
modelIrrev = changeRxnBounds(modelIrrev,'EnzCon_C_Rxn',0.08,'l' ); % tuned by Trinankur
EnzCon_MetInd=find(ismember(modelIrrev.mets,'EnzCon_Met[c]'));
Rxn_WithEnzConInd=find(ismember(modelIrrev.rxns,kcat_rxns_irrev));
selExc=findExcRxns(modelIrrev);
ExRxnInd=find(selExc);
BiomassRxnInd=find(ismember(modelIrrev.rxns,biomass));
Solutions=zeros(length(modelIrrev.rxns),sample_size);
EnzConSamples=zeros(length(modelIrrev.rxns)-1,sample_size);
fluxSum = zeros(length(modelIrrev.mets),sample_size);           %Added by VIKASH & KANNAN
RandomEnzConCoeff=zeros(length(modelIrrev.rxns),sample_size);
sort_kcat_mw = sortrows(kcat_mw);
for i=1:1:sample_size
    RandomEnzConCoeff(:,i)=randsample(kcat_mw,length(modelIrrev.rxns),true);
%     RandomEnzConCoeff(:,i)=randsample(sort_kcat_mw(1:max(find(sort_kcat_mw == prctile(sort_kcat_mw,90)))),length(modelIrrev.rxns),true);
%     RandomEnzConCoeff(:,i)=median(kcat_mw);
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
%         if(isempty(modelIrrev1.S(EnzCon_MetInd,Rxn_WithEnzConInd(k))))==0
            modelIrrev1.S(EnzCon_MetInd,Rxn_WithEnzConInd(k))=kcat_mw1(find(ismember(kcat_rxns_irrev,modelIrrev1.rxns(Rxn_WithEnzConInd(k)))));
%         end
    end
    %% Removing EnzCon coeff for exchange reactions
    modelIrrev1.S(EnzCon_MetInd,ExRxnInd)=0;
    modelIrrev1.S(EnzCon_MetInd,find(ismember(modelIrrev1.rxns,'EnzCon_C_Rxn')))=-1;
%     changeCobraSolver('gurobi7');
    solution=optimizeCbModel(modelIrrev1,'max','one');
%     solution=optimizeCbModel(modelIrrev1,'max');
    if(isempty(solution.x))==0
        if rem(j,500) == 0
            fprintf('%d Solutions Generated\n\n',j)
        end
        Solutions(:,j)=solution.x;
        fs = sum(abs((modelIrrev1.S)*(diag(solution.x)))' * 0.5)';           %Added by VIKASH & KANNAN
        fluxSum(:,j)=fs;
%         ShadowPrices(:,j)=solution.y;
%         ReducedCosts(:,j)=solution.w;
%         EnzConSamples(:,j)=modelIrrev1.S(EnzCon_MetInd,1:end-1)';
    end
end
valid_sol = find(Solutions(BiomassRxnInd,:)>0);
Solutions = Solutions(:,valid_sol);
fluxSum = fluxSum(:,valid_sol);
% EnzConSamples = EnzConSamples(1:end,valid_sol);
% ShadowPrices = ShadowPrices(1:end,valid_sol);
% ReducedCosts = ReducedCosts(1:end,valid_sol);
end
function [model,rxnIDexists] = addReaction_local(model,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate)
%addReaction Add a reaction to the model or modify an existing reaction
%
% model = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate)
% model = addReaction(model,rxnName,rxnFormula)
%
%INPUTS
% model             COBRA model structure
% rxnName           Reaction name abbreviation (i.e. 'ACALD')
%                   (Note: can also be a cell array {'abbr','name'}
% metaboliteList    Cell array of metabolite names or alternatively the
%                   reaction formula for the reaction
% stoichCoeffList   List of stoichiometric coefficients (reactants -ve,
%                   products +ve), empty if reaction formula is provided
%
%OPTIONAL INPUTS
% revFlag           Reversibility flag (Default = true)
% lowerBound        Lower bound (Default = 0 or -vMax)
% upperBound        Upper bound (Default = vMax)
% objCoeff          Objective coefficient (Default = 0)
% subSystem         Subsystem (Default = '')
% grRule            Gene-reaction rule in boolean format (and/or allowed)
%                   (Default = '');
% geneNameList      List of gene names (used only for translation from 
%                   common gene names to systematic gene names)
% systNameList      List of systematic names
% checkDuplicate    Check S matrix too see if a duplicate reaction is
%                   already in the model (Deafult true)
%
%OUTPUTS
% model             COBRA model structure with new reaction
% rxnIDexists       Empty if the reaction did not exist previously, or if
%                   checkDuplicate is false. Otherwise it contains the ID
%                   of an identical reaction already present in the model.
%
% Examples:
%
% 1) Add a new irreversible reaction using the formula approach
%
%    model = addReaction(model,'newRxn1','A -> B + 2 C')
%
% 2) Add a the same reaction using the list approach
%
%    model = addReaction(model,'newRxn1',{'A','B','C'},[-1 1 2],false);
%
% Markus Herrgard 1/12/07
%
% Modified the check to see if duplicate reaction already is in model by
% using S matrix coefficients to be able to handle larger matricies
% Richard Que 11/13/2008


parseFormulaFlag = false;
rxnIDexists = [];

if iscell(rxnName)&&length(rxnName)>1
    rxnNameFull = rxnName{2};
    rxnName = rxnName{1};
end

% Figure out if reaction already exists
nRxns = length(model.rxns);
if (sum(strcmp(rxnName,model.rxns)) > 0)
    warning('Reaction with the same name already exists in the model');
    [tmp,rxnID] = ismember(rxnName,model.rxns);
    oldRxnFlag = true;
else
    rxnID = nRxns+1;
    oldRxnFlag = false;
end

% Figure out what input format is used
if (nargin < 4)
    if (~iscell(metaboliteList))
        parseFormulaFlag = true;
    else
        error('Missing stoichiometry information');
    end
else
    if isempty(stoichCoeffList)
        parseFormulaFlag = true;
    else
        if (length(metaboliteList) ~= length(stoichCoeffList))
            error('Incorrect number of stoichiometric coefficients provided');
        end
    end
end

% Reversibility
if (nargin < 5 | isempty(revFlag))
    if (oldRxnFlag)
        revFlag = model.rev(rxnID);
    else
        revFlag = true;
    end
end

% Parse formula
if (parseFormulaFlag)
    rxnFormula = metaboliteList;
    [metaboliteList,stoichCoeffList,revFlag] = parseRxnFormula(rxnFormula);
end

% Missing arguments
if (nargin < 6 | isempty(lowerBound))
    if (oldRxnFlag)
        lowerBound = model.lb(rxnID);
    else
        if (revFlag)
            lowerBound = min(model.lb);
            if isempty(lowerBound)
                lowerBound=-1000;
            end
        else
            lowerBound = 0;
        end
    end
end
if (nargin < 7 | isempty(upperBound))
    if (oldRxnFlag)
        upperBound = model.ub(rxnID);
    else
        upperBound = max(model.ub);
        if isempty(upperBound)
            upperBound=1000;
        end
    end
end
if (nargin < 8 | isempty(objCoeff))
    if (oldRxnFlag)
        objCoeff = model.c(rxnID);
    else
        objCoeff = 0;
    end
end
if (nargin < 9)
    if (oldRxnFlag) && (isfield(model,'subSystems'))
        subSystem = model.subSystems{rxnID};
    else
        subSystem = '';
    end
end
if  isempty(subSystem)
    if (oldRxnFlag) && (isfield(model,'subSystems'))
        subSystem = model.subSystems{rxnID};
    else
        subSystem = '';
    end
end
if (nargin < 10) && (isfield(model,'grRules'))
    if (oldRxnFlag)
        grRule = model.grRules{rxnID};
    else
        grRule = '';
    end
end

if (~exist('checkDuplicate','var'))
   checkDuplicate=true;
end

nMets = length(model.mets);
Scolumn = sparse(nMets,1);

modelOrig = model;

% Update model fields
model.rxns{rxnID,1} = rxnName;
if (revFlag)
    model.rev(rxnID,1) = 1;
else
    model.rev(rxnID,1) = 0;
end
model.lb(rxnID,1) = lowerBound;
model.ub(rxnID,1) = upperBound;
model.c(rxnID,1) = objCoeff;

if (isfield(model,'rxnNames'))
    if exist('rxnNameFull','var')
        model.rxnNames{rxnID,1} = rxnNameFull;
    else
        model.rxnNames{rxnID,1} = model.rxns{rxnID};
    end
end
if (isfield(model,'subSystems'))
    model.subSystems{rxnID,1} = subSystem;
end
if isfield(model,'rxnNotes')
    model.rxnNotes{rxnID,1} = '';
end
if isfield(model,'confidenceScores')
    model.confidenceScores{rxnID,1} = '';
end
if isfield(model,'rxnReferences')
    model.rxnReferences{rxnID,1} = '';
end
if isfield(model,'rxnECNumbers')
    model.rxnECNumbers{rxnID,1} = '';
end


% Figure out which metabolites are already in the model
[isInModel,metID] = ismember(metaboliteList,model.mets);

nNewMets = sum(~isInModel);

% Construct S-matrix column
newMetsCoefs=zeros(0);
for i = 1:length(metaboliteList)
    if (isInModel(i))
        Scolumn(metID(i),1) = stoichCoeffList(i);
    else
        warning(['Metabolite ' metaboliteList{i} ' not in model - added to the model']);
        Scolumn(end+1,1) = stoichCoeffList(i);
        model.mets{end+1,1} = metaboliteList{i};
        newMetsCoefs(end+1) = stoichCoeffList(i);
        if (isfield(model,'metNames'))      %Prompts to add missing info if desired
            model.metNames{end+1,1} = regexprep(metaboliteList{i},'(\[.+\]) | (\(.+\))','') ;
            warning(['Metabolite name for ' metaboliteList{i} ' set to ' model.metNames{end}]);
%             model.metNames(end) = cellstr(input('Enter complete metabolite name, if available:', 's'));
        end
        if (isfield(model,'metFormulas'))
            model.metFormulas{end+1,1} = '';
            warning(['Metabolite formula for ' metaboliteList{i} ' set to ''''']);
%             model.metFormulas(end) = cellstr(input('Enter metabolite chemical formula, if available:', 's'));
        end
        if isfield(model,'metChEBIID')
            model.metChEBIID{end+1,1} = '';
        end
        if isfield(model,'metKEGGID')
            model.metKEGGID{end+1,1} = '';
        end
        if isfield(model,'metPubChemID')
            model.metPubChemID{end+1,1} = '';
        end
        if isfield(model,'metInChIString')
            model.metInChIString{end+1,1} = '';
        end
        if isfield(model,'metCharge')
            model.metCharge(end+1,1) = 0;
        end
    end
end

%printLabeledData(model.mets,Scolumn,1);

if isfield(model,'b')
    model.b = [model.b;zeros(length(model.mets)-length(model.b),1)];
end

% if ~oldRxnFlag, model.rxnGeneMat(rxnID,:)=0; end

if (isfield(model,'genes'))
    if (nargin < 11)
        model = changeGeneAssociation(model,rxnName,grRule);
    else
        model = changeGeneAssociation(model,rxnName,grRule,geneNameList,systNameList);
    end
end

% Figure out if the new reaction already exists
rxnInModel=false;
if (nNewMets > 0) && isempty(find(newMetsCoefs == 0, 1))
    Stmp = [model.S;sparse(nNewMets,nRxns)];
else
    Stmp = model.S;
    if (checkDuplicate)
       if size(Stmp,2)<6000
           tmpSel = all(repmat((Scolumn),1,size(Stmp,2)) == (Stmp));
           rxnIDexists = full(find(tmpSel));
           if (~isempty(rxnIDexists))
               rxnIDexists=rxnIDexists(1);
               rxnInModel = true;
           end
       else
           for i=1:size(Stmp,2)
               if(Scolumn==Stmp(:,i))
                   rxnInModel=true;
                   rxnIDexists=i;
                   break
               end
           end
       end
    end
end

if (rxnInModel)
    warning(['Model already has the same reaction you tried to add: ' modelOrig.rxns{rxnIDexists}]);
    model = modelOrig;
else
    if (oldRxnFlag)
        model.S = Stmp;
        model.S(:,rxnID) = Scolumn;
    else
        model.S = [Stmp Scolumn];
    end
%     printRxnFormula(model,rxnName);
end
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