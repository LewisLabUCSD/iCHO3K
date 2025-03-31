function Nsnp = fast_snp(S,lb,ub,vTol)
% Fast Sparse Nullspace Pursuit
%
% Solves a series of LP problems of the form
%
%    min      ~& ||v||_1 \\
%    s.t.     ~& S v = 0 \\
%             ~& x^T v \geq epsilon \\
%             ~& or x^T v \leq -epsilon \\
%             ~& lb \leq v \leq ub \\
%
% USAGE:
%              Nsnp = fastSNP(S, lb, ub, vTol)
%
% INPUTS:
%              S:  m x n Stoichiometric matrix
%              lb: n x 1 Lower bounds
%              ub: n x 1 Upper bounds
%
% OPTIONAL INPUTS:
%              vTol: zero flux tolerance (double)
%
% OUTPUT:
%              Nsnp: n x (n-rank(S)) matrix with a sparse null space basis
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Define solver parameters
if nargin<3
    disp('Not enough input parameters supported'); 
    Nsnp = [];
    return;
elseif nargin<4
    vTol = 1e-6;
end

% Initialization of SNP parameters
nullity = size(null(full(S)),2);
Nsnp    = [];
P_N     = 0;
epsilon = 1e-1;

% Main loop iterates until the desired number of basis vectors is reached
% or, until there is no feasible solution left
for jx = 1:nullity
    
    %Initialize non-negative random weights
    weights = rand(1,size(S,2));
    vopt    = [];
    P_NT    = eye(size(S,2))-P_N'*P_N;
    
    % Build positive SNV problem considering the fast-SNP mode
    x          = (1+weights)*P_NT;
    SNVproblem = buildSNVproblem(S,lb,ub,x,epsilon);
    solPos     = solveCobraLP(SNVproblem);
    
    % Build negative SNV problem considering the fast-SNP mode
    SNVproblem = buildSNVproblem(S,lb,ub,-x,epsilon);
    solNeg     = solveCobraLP(SNVproblem);
    
    % Get SNV solutions
    if (solPos.stat==1)
        vopt = solPos.full(1:size(S,2));
        vopt(abs(vopt)<vTol) = 0;
        posLength = sum(vopt~=0);
    else
        posLength = Inf;
    end
    if (solNeg.stat==1)
        vopt = [vopt,solNeg.full(1:size(S,2))];
        vopt(abs(vopt(:,end))<vTol,end) = 0;
        negLength = sum(vopt(:,end)~=0);
    else
        negLength = Inf;
    end

    % Choose shortest (valid) solution
    if ~isempty(vopt)
        if ~isinf(posLength) && ~isinf(negLength)
            [~,ix] = min([posLength,negLength]);
            Nsnp = [Nsnp,vopt(:,ix)];            % Choose shortest solution
        else
            Nsnp = [Nsnp,vopt];
        end
    else
        break;
    end

    % Re-build projetion onto col(x)
    P_N = orth(Nsnp)';
end
