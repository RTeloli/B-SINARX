function [y_pre, err] = closedloopVarDelay3(mathFunctions, A, delays, arecur)
% Closed-loop simulation with a different delay for each variable.
%   A      : data (time x variables); transposed automatically if needed
%   delays : delay for each variable (length = number of columns of A)
%   arecur : logical mask of the recurrent (fed-back) variables

if size(A, 1) < size(A, 2), A = A.'; end
leng   = size(A, 1);
arecur = logical(arecur);
maxd   = max(delays);

% keep(r, j) is true for the last delays(j) rows of column j of the window
keep  = (1:maxd)' > maxd - delays(:)';
scale = max(abs(A), [], "all");

% Initialise with the true values up to the maximum delay
y_pre = zeros(leng, nnz(arecur));
y_pre(1:maxd, :) = A(1:maxd, arecur);
err   = zeros(leng, nnz(arecur));

for t = maxd+1 : leng
    % Past maxd steps: true exogenous inputs, predicted recurrent states
    W = A(t-maxd : t-1, :);
    W(:, arecur) = y_pre(t-maxd : t-1, :);

    y = predictFunction(mathFunctions, W(keep).');
    y_pre(t, :) = y;
    err(t, :)   = A(t, arecur) - y;

    % Divergence checks
    if any(isnan(y))
        y_pre(t:end, :) = NaN;
        fprintf("NaN in prediction at step %d: divergent solution!\n", t);
        return
    elseif max(abs(y)) / scale > 20
        y_pre(t:end, :) = NaN;
        fprintf("Prediction exceeds 20x data range at step %d: divergent solution!\n", t);
        return
    end
end
end