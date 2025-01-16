function y = sampleSet(n, k, varargin)
%RANDSAMPLE Random sample, with or without replacement.
%   Y = RANDSAMPLE(N,K) returns Y as a column vector of K values sampled
%   uniformly at random, without replacement, from the integers 1:N.
%
%   Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at random,
%   without replacement, from the values in the vector POPULATION.  Y is a
%   vector of the same type as POPULATION.
% 
%   NOTE:  When POPULATION is a numeric vector containing only non-negative
%   integer values, and it might have length 1, use either
%
%      Y = DATASAMPLE(POPULATION,K,'Replace',false)
%   or
%      Y = POPULATION(RANDSAMPLE(LENGTH(POPULATION),K))
%

nargs = nargin;

% Process the stream argument, if present
defaultStream = isnumeric(n) || ~iIsRNG(n); % simple test first for speed
if defaultStream
    s = RandStream.getGlobalStream();
else
    % shift right to drop s from the argument list
    nargs = nargs - 1;
    s = n;
    n = k;
    if nargs > 1
        k = varargin{1};
        varargin(1) = [];
    end
end

if isscalar(n) && isnumeric(n) && (round(n) == n) && (n >= 0)
    havePopulation = false;
    population = [];
else
    havePopulation = true;
    population = n;
    n = numel(population);
    if ~isvector(population)
        error(message('stats:randsample:BadPopulation'));
    end
end

if nargs < 3
    replace = false;
else
    replace = varargin{1};
end

if nargs < 4
    w = [];
else
    w = varargin{2};
    if ~isempty(w)
        if length(w) ~= n
            error(message('stats:randsample:InputSizeMismatch', n));
        else
            sumw = sum(w);
            if ~(sumw > 0) || ~all(w>=0) % catches NaNs
                error(message('stats:randsample:InvalidWeights'));
            end
            p = w(:)' / sumw;
        end
    end
end

switch replace
    
    % Sample with replacement
    case {true, 'true', 1}
        if n == 0
            if k == 0
                y = zeros(0,1);
            else
                error(message('stats:randsample:EmptyPopulation'));
            end
        elseif isempty(w)
            y = randi(s,n,k,1);
            
        else
            edges = min([0 cumsum(p)],1); % protect against accumulated round-off
            edges(end) = 1; % get the upper edge exact
            [~, ~, y] = histcounts(rand(s,k,1),edges);
            % Type of Y should match edges
            y = cast(y, "like", edges);
        end
        
        % Sample without replacement
    case {false, 'false', 0}
        if k > n
            error(message('stats:randsample:SampleTooLarge', n));
        end
        
        if isempty(w)
            % If the sample is a sizable fraction of the population,
            % just randomize the whole population (which involves a full
            % sort of n random values), and take the first k.
            if 4*k > n
                rp = randperm(s,n);
                y = rp(1:k);
                
                % If the sample is a small fraction of the population, a full sort
                % is wasteful.  Repeatedly sample with replacement until there are
                % k unique values.
            else
                x = false(1,n); % flags
                sumx = 0;
                while sumx < floor(k) % prevent infinite loop when 0<k<1
                    x(randi(s,n,1,k-sumx)) = true; % sample w/replacement
                    sumx = nnz(x); % count how many unique elements so far
                end
                y = find(x);
                y = y(randperm(s,k));
            end
        else
            error(message('stats:randsample:NoWeighting'));
        end
    otherwise
        error(message('stats:randsample:BadReplaceValue'));
end

if havePopulation
    y = population(y);
else
    y = y(:);
end