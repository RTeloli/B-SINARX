function xi = trimmedSR3(theta, dX, k)

lambda_rel = logspace(-3, 1, 100);
lambdan = norm(dX);
if   abs(lambdan)<eps
    lambdan = eps;
end
counter = 1; 
while true
    for i = 1:length(lambda_rel)
        lambda = lambda_rel(i)*lambdan;
        xi  = reduceXi(theta, dX', lambda);
        qqq = nnz(xi);
        counter = counter+1;
        if qqq<=k
            break 
        end
    end
    if qqq==k || counter>10000 || abs(lambdan)<eps 
        break
    elseif qqq>=k
        lambdan = 1.1*lambdan;
    elseif qqq<=k
        lambdan = 0.9*lambdan;
    end
end

end
function xi  = reduceXi(theta_x, y, lambda)
[n, m] = size(theta_x);
if n<m
    theta_x=theta_x';
end
[n, m] = size(y);
if n<m
    y=y';
end

n_basis = size(theta_x,2);
[~, n_outputs] = size(y) ;
xi = zeros(n_basis, n_outputs);
for i=1:n_outputs
    [~, w] = sr3(theta_x,y(:, i),  "lam", lambda(i));
    xi(:, i) = w;
end
end