function obj = trainSTLS(obj, theta_norm, dX)
obj.trainedalgorithm = "STLS";

trainRatio = obj.trainRatio;
[theta_train,dX_train,theta_test, dX_test ] = splitData(theta_norm, dX, trainRatio) ;

lambdan = vecnorm(dX, 2);
if   any(abs(lambdan)<eps)
    lambdan = eps+lambdan;
end

lambda_0            = zeros(obj.mode.n, 1);
xi_0                = reduceXi( theta_train,dX_train,lambda_0);

[theta_train_reduced, theta_test_reduced]  = reduceTheta(theta_train, theta_test, xi_0);

n =500;
lambda_rel = logspace(-7, 2, n );

lmb = lambda_rel'*lambdan;

spar = zeros(n, size(lmb, 2));
nnormxi = zeros(n, size(lmb, 2));
clc
fprintf("Running Optimization\n")

parfor i = 1:n
    xi{i}         = reduceXi(theta_train_reduced, dX_train, lmb(i, :));
    spar(i, :)    = size(xi{i}, 1) - sum(xi{i}==0);
    nnormxi(i, :) = vecnorm(xi{i}, 1);

end
[~, kk, loss_1] = RICc(theta_test_reduced, dX_test , xi);


for i = 1:obj.mode.n
    use = spar(:, i)>1;

    f1         = loss_1(use, i) /max(loss_1(use, i));
    f2         = spar(use, i) ;%/max(spar(:, i));

    if any(isnan(f1)) || any(isnan(f2))
        maxIndex = 1;
        f1  = loss_1(:, i);
        f2  = spar(:, i);
    elseif cov(f1)<1e-4

        maxIndex = length(f2);

    else

        maxIndex   =  traingleOPt(f1, f2 );
    end
    % f1 = 1-f1;
    if obj.pltOption
        paretoPlot(f1, f2, maxIndex, ...
            "title",  "STLS", ...
            "save", obj.savePlots, ...
            "var", obj.outputVars{i}, ...
            "noise", obj.noise)
    end
    lmb_try(i) = lmb(maxIndex, i);
end

xi_c = reduceXi(theta_train, dX_train,lmb_try);

obj.Xi = xi_c;

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
% xi = (theta_x'*theta_x + 0.1*eye)\theta_x'*y;
[~, n_outputs] = size(y) ;
xi = zeros(n_basis, n_outputs);
for i=1:n_outputs
    xi(:, i) = sparsifyDynamics(theta_x,y(:, i),  lambda(i), 1);
end

end

