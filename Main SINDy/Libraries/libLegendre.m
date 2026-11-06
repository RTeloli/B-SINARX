function [sol, vars, mainString, l] = libLegendre(x, p, vars, q)
% LIBLEGENDRE  Build an orthogonal (Legendre) feature library up to order p,
%              including multivariate cross terms via hyperbolic truncation.
%
%   x    : matrix with one dimension = length(vars) (variables),
%          the other = number of samples. Orientation is auto-detected.
%          Each column of x (after orientation fix) should be normalized to [-1, 1].
%   p    : max total Legendre order (sum of per-variable orders)
%   vars : cell array of variable name strings
%   q    : hyperbolic truncation parameter, 0 < q <= 1
%          q = 1  -> full tensor product up to total order p (all combinations)
%          q < 1  -> sparser set, penalizes high-order cross terms
%          (default: 1, if not provided)
%
%   sol        : n x M library matrix
%   mainString : M x 1 string array of term labels (e.g. "legN(1,x1)*legN(2,u1)")
%   l          : M x nvars matrix, l(k,j) = Legendre order of vars{j} in term k

if nargin < 4
    q = 1;
end

nvars = length(vars);

% --- Auto-detect orientation ---
[r, c] = size(x);
if c == nvars
    % already n x nvars
elseif r == nvars
    x = x';
else
    error('libLegendre:sizeMismatch', ...
        'Neither dimension of x (%dx%d) matches length(vars)=%d.', r, c, nvars);
end

n = size(x, 1);

% --- Build hyperbolic-truncated multi-index set (each row = orders per variable) ---
l = hyperbolicIndexSet(nvars, p, q);   % M x nvars, includes the all-zero (constant) row
M = size(l, 1);

% --- Precompute numeric orthonormal Legendre basis per variable, orders 0..p ---
Puni = cell(nvars, 1);
for j = 1:nvars
    Puni{j} = legendreRecurrence(x(:,j), p);   % n x (p+1)
end

% --- Assemble library by multiplying univariate bases per multi-index row ---
sol = ones(n, M);
mainString = strings(M, 1);
for k = 1:M
    labelParts = strings(0);
    for j = 1:nvars
        ord = l(k, j);
        if ord > 0
            sol(:,k) = sol(:,k) .* Puni{j}(:, ord+1);
            labelParts(end+1) = sprintf("legN(%d,%s)", ord, vars{j}); %#ok<AGROW>
        end
    end
    if isempty(labelParts)
        mainString(k) = "1";
    else
        mainString(k) = strjoin(labelParts, "*");
    end
end

end

function P = legendreRecurrence(x, maxOrder)
% Numeric orthonormal Legendre basis, N x (maxOrder+1), orders 0..maxOrder
    N = length(x);
    P = zeros(N, maxOrder+1);
    P(:,1) = 1;
    if maxOrder >= 1
        P(:,2) = x;
    end
    for n = 2:maxOrder
        P(:,n+1) = ((2*n-1)*x.*P(:,n) - (n-1)*P(:,n-1)) / n;
    end
    for n = 0:maxOrder
        P(:,n+1) = P(:,n+1) * sqrt((2*n+1)/2);
    end
end

function idxSet = hyperbolicIndexSet(d, p, q)
% Generate multi-indices (rows) of length d with each entry in [0,p],
% total order <= p, filtered by hyperbolic (q-norm) truncation:
%   (sum(alpha.^q))^(1/q) <= p
% Includes the all-zero row (constant term) first.
    idxSet = [];
    % Generate all combinations with total order <= p via recursive expansion
    idxSet = genIndices(d, p);
    % Apply hyperbolic truncation
    normVals = sum(idxSet.^q, 2).^(1/q);
    keep = normVals <= p + 1e-9;
    idxSet = idxSet(keep, :);
    % Ensure all-zero row is first
    isZeroRow = all(idxSet == 0, 2);
    idxSet = [idxSet(isZeroRow,:); idxSet(~isZeroRow,:)];
end

function idxSet = genIndices(d, p)
% All d-tuples of nonnegative integers with sum <= p
    if d == 1
        idxSet = (0:p)';
        return;
    end
    idxSet = [];
    for total = 0:p
        idxSet = [idxSet; genIndicesFixedSum(d, total)]; %#ok<AGROW>
    end
end

function idxSet = genIndicesFixedSum(d, s)
% All d-tuples of nonnegative integers summing exactly to s
    if d == 1
        idxSet = s;
        return;
    end
    idxSet = [];
    for first = 0:s
        rest = genIndicesFixedSum(d-1, s-first);
        idxSet = [idxSet; first*ones(size(rest,1),1), rest]; %#ok<AGROW>
    end
end