function E_C = calc_conn_evidence(model, E_X)

% S matrix is binarized to indicate metabolite participation in each reaction
Sbin = double(model.S ~= 0);

% Adjacency matrix (i.e., binary reaction connectivity); the connection between
% a reaction and itself is ignored by subtracting the identity matrix
A = full(double(Sbin' * Sbin ~= 0));
A = A - eye(size(A));

% Influence matrix; describes the divided connectivity of reactions --
% e.g., if R1 is connected to 4 reactions, its influence on each of those
% reactions would be 0.25
I = A ./ repmat(sum(A, 2), 1, size(A, 2));

% Weighted influence matrix; each reaction's influence on others is
% weighted by its expression score
WI = repmat(E_X, 1, size(A, 2)) .* I;

% Connectivity score; sum the influence of all connected reactions
E_C = sum(WI)';