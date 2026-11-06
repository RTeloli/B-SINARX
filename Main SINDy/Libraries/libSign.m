function [sol, vars, mainString, l] = libSign(x, p, vars)
% LIBSIGN  Build a sign-based (Coulomb friction type) feature library
%          up to power p: sign(x), sign(x)*x, sign(x)*x^2, ..., sign(x)*x^p
%          per variable. Linear-in-parameters (coefficients enter linearly).
%
%   x    : matrix with one dimension = length(vars) (variables),
%          the other = number of samples. Orientation is auto-detected.
%   p    : max power multiplying sign(x) (p=0 gives only sign(x) terms)
%   vars : cell array of variable name strings
%
%   sol        : n x (1+p*nvars) library matrix (col 1 = constant term "1")
%   mainString : (1+p*nvars) x 1 string array of term labels
%   l          : (1+p*nvars) x nvars matrix, l(k,j) = power of vars{j} in
%                term k (sign() itself is not tracked in l, only the power)

nvars = length(vars);

% --- Auto-detect orientation: rows should be samples, cols = variables ---
[r, c] = size(x);
if c == nvars
    % already n x nvars, do nothing
elseif r == nvars
    x = x';   % transpose to n x nvars
else
    error('libSign:sizeMismatch', ...
        'Neither dimension of x (%dx%d) matches length(vars)=%d.', r, c, nvars);
end

n = size(x, 1);
nTerms = 1 + p*nvars;

% --- Build library matrix and exponent matrix ---
sol = ones(n, 1);
l   = zeros(nTerms, nvars);   % row per term, col per variable
rowIdx = 2;  % row 1 (constant term) stays all-zero

for i = 0:p-1
    block = zeros(n, nvars);
    for j = 1:nvars
        if i == 0
            block(:,j) = sign(x(:,j));              % sign(x)
        else
            block(:,j) = sign(x(:,j)) .* x(:,j).^i;  % sign(x)*x^i
        end
        l(rowIdx, j) = i;   % track the power (sign itself implied by this library)
        rowIdx = rowIdx + 1;
    end
    sol = [sol, block];
end

% --- Build matching string labels ---
mainString = "1";
for i = 0:p-1
    for j = 1:nvars
        if i == 0
            str = sprintf("sign(%s)", vars{j});
        else
            sup = toSuperscript(num2str(i));
            str = sprintf("sign(%s)*%s%s", vars{j}, vars{j}, sup);
        end
        mainString = [mainString, str];
    end
end
mainString = mainString';

end

function out = toSuperscript(str)
map = containers.Map(...
    {'0','1','2','3','4','5','6','7','8','9'}, ...
    {'⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹'});
out = '';
for k = 1:length(str)
    if isKey(map, str(k))
        out = [out map(str(k))];
    else
        out = [out str(k)];
    end
end
end