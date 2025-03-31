function A = normMatrixEntries(A)
% Normalized (sparse) matrix
%
% Normalizes each column by the largest entry
%
% USAGE:
%              A = normMatrixEntries(A)
%
% INPUTS:
%              A - m x n  sparse basis matrix
%
% OUTPUT:
%              A - m x n  normalized sparse basis matrix
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

A(~A) = inf;                                   % Change zeros for inf
A     = bsxfun(@rdivide,A,min(abs(A)));        % Extract topological features by diving for smallest non-zero flux
A(~isfinite(A)) = 0;                           % Transform inf entries to the original zeros
