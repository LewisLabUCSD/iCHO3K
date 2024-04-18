function [A,index] = linearization_index(A,flag)

% USAGE:
% % [A,index] = linearization_index(A,flag)

% INPUTS:
% % A:          a cell array of cell arrays
% % flag:       'cols', if each element is of the form {nx1 cell}
                % % 'rows' if each element is of the form {1xn cell}

% OUTPUTS:
% % A:          a cell array where subcell arrays have been opened and concatenated
% % index:      a numeric describing index mapping from old A to new A

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

cnt = 0;
index = zeros(1,1);
if iscell(A{1,1})
    for i=1:length(A)
        k = length(A{i,1});
        for j=1:k
            cnt = cnt+1;
            index(cnt,1) = i;
            if strcmp(flag,'rows')
                B{cnt,1} = A{i,1}{1,j};
            else
                B{cnt,1} = A{i,1}{j,1};
            end
        end
    end
else
    for i=1:length(A)
        k = length(A{i,1});
        for j=1:k
            cnt = cnt+1;
            index(cnt,1) = i;
            if strcmp(flag,'rows')
                B(cnt,1) = A{i,1}(1,j);
            else
                B(cnt,1) = A{i,1}(j,1);
            end
        end
    end
end
A = B;