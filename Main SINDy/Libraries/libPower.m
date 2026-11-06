function [sol, vars, mainString, l] = libPower(x, p, vars)
% LIBPOWER  Build a polynomial feature library up to power p
%
%   x    : matrix with one dimension = length(vars) (variables),
%          the other = number of samples. Orientation is auto-detected.
%   p    : max polynomial power
%   vars : cell array of variable name strings
%
%   sol        : n x (1+p*nvars) library matrix
%   mainString : (1+p*nvars) x 1 string array of term labels
%   l          : (1+p*nvars) x nvars matrix, l(k,j) = power of vars{j} in term k

nvars = length(vars);

% --- Auto-detect orientation: rows should be samples, cols = variables ---
[r, c] = size(x);
if c == nvars
    % already n x nvars, do nothing
elseif r == nvars
    x = x';   % transpose to n x nvars
else
    error('libPower:sizeMismatch', ...
        'Neither dimension of x (%dx%d) matches length(vars)=%d.', r, c, nvars);
end

n = size(x, 1);
nTerms = 1 + p*nvars;

% --- Build library matrix and exponent matrix ---
sol = ones(n, 1);
l   = zeros(nTerms, nvars);   % row per term, col per variable

rowIdx = 2;  % row 1 (constant term) stays all-zero
for i = 1:p
    sol = [sol, x.^i];
    for j = 1:nvars
        l(rowIdx, j) = i;
        rowIdx = rowIdx + 1;
    end
end

% --- Build matching string labels ---
mainString = "1";
for i = 1:p
    for j = 1:nvars
        if i ~= 1
            sup = toSuperscript(num2str(i));
            str = strcat(vars{j}, sup);
        else
            str = vars{j};
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