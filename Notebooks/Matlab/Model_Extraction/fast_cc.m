function [inactive] = fast_cc(gurobi_model,rxn_indices)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    rxn_indices = 1:numel(gurobi_model.lb);
    rxn_indices = rxn_indices';
end

tol = 1e-8;

inactive = false(size(gurobi_model.lb));

%done = false;
gurobi_model.obj(rxn_indices) = 1;
gurobi_model.modelsense = 'max';
params.outputflag = 0;
params.FeasibilityTol = 1e-9;
params.OptimalityTol = 1e-6;
%nx0 = sum(gurobi_model.obj);
r1 = rxn_indices;
while ~isempty(r1)
    gurobi_model.obj(1:end) = 0;
    gurobi_model.obj(r1(1)) = 1;
    r = gurobi(gurobi_model,params);
    x = abs(r.x);
    active = find(x >= tol);
    r1(ismember(r1,active)) = [];
    
    if abs(r.objval) <= tol
        gurobi_model.modelsense = 'min';
        r = gurobi(gurobi_model,params);
        x = abs(r.x);
        active = find(x >= tol);
        r1(ismember(r1,active)) = [];
        if abs(r.objval) <= tol
            inactive(logical(gurobi_model.obj)) = true;
            r1(ismember(r1,find(gurobi_model.obj))) = [];
        end
        gurobi_model.modelsense = 'max';
    end
end
inactive = inactive(rxn_indices);

end

