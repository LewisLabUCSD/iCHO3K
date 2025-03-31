function isFeasible = looplessCheck(points,sample)
% Checks the loopless condition topologically
%
% Checks loopless condition of a series of points for a given sample structure
%
% USAGE:
%              isFeasible = checkFeasibility(points, sample)
%
% INPUTS:
%              points:  n x k Matrix of points
%              sample (structure):  (the following fields are required - others can be supplied)
%                                    * loopMatrix - `ni x p` Loop law matrix
%                                    * nnzEntries - `1 x p` Matrix of non-zero entries of the loop law matrix
%                                    * intRxns - `ni x 1` Matrix with indexes of internal rxns
%
% OUTPUT:
%              isFeasible:  1 x k boolean vector
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

loopyFlux = false(1,size(points,2));
for ix = 1:numel(sample.nnzEntries)
    loopyFlux = (loopyFlux | bsxfun(@eq,abs(sum(bsxfun(@times,sign(points(sample.intRxns,:)),sample.loopMatrix(:,ix)))),sample.nnzEntries(ix)));
end
isFeasible = ~loopyFlux;
