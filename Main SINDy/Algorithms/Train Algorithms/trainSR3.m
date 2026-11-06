function obj = trainSR3(obj, theta_norm, dX)
obj.trainedalgorithm = "SR3";
lambdan = vecnorm(dX, 2);
if   any(abs(lambdan)<eps)
    lambdan = eps+lambdan;
end
trainRatio = obj.trainRatio;
[theta_train,dX_train,theta_test, dX_test ] = splitData(theta_norm, dX, trainRatio) ;

lambda_0            = eps*ones(obj.mode.n, 1);
xi_0                = reduceXi( theta_train,dX_train,lambda_0);

[theta_train_reduced, theta_test_reduced]  = reduceTheta(theta_train, theta_test, xi_0);

% xi_0_cpy            = xi_0;
% [~,inx]             = maxk(xi_0, 1);

% for i = 1:obj.mode.n
%     xi_0_cpy(inx(:, i), i) =NaN;
% end
% xi_0_cpy(xi_0_cpy ==0) = NaN;
n = 5000; 
lambda_rel = logspace(-6, 3, n );

lmb = lambda_rel'*lambdan;
% lmb =   lmbmin' + (lmbmax -lmbmin )'*linspace(0, 1, n );
% lmb = eps + lmb';

spar = zeros(n, size(lmb, 2));
nnormxi = zeros(n, size(lmb, 2));
loss_1 = zeros(n, size(lmb, 2));

clc
fprintf("Running Optimization\n")

parfor i = 1:n
    xi{i}      = reduceXi(theta_train_reduced, dX_train, lmb(i, :));
    spar(i, :) = size(xi{i}, 1) - sum(xi{i}==0);
    nnormxi(i, :) = vecnorm(xi{i}, 1); 
    y_pred_1       = theta_test_reduced*xi{i};
     loss_1(i,:) = 0.5*vecnorm(dX_test - y_pred_1, 2);
end
[~, kk, loss] = RICc(theta_test_reduced, dX_test , xi);


for i = 1:obj.mode.n
    if max(loss(:, i)) == 0
            maxIndex = 1;
            f1  = loss(:, i);
            f2  = spar(:, i);
    else
        f1        =  loss(:, i) /max(loss(:, i));
        f2        =  spar(:, i);%/max(spar(:, i));
        % f1 = 1-f1;
        maxIndex   =  traingleOPt(f1, f2 );
    end
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
obj =  obj.makeEquation();
% obj.MathFuncts = generateMathFunction(obj);

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
    [~, w] = sr3(theta_x,y(:, i),  "lam", lambda(i));
    xi(:, i) = w; 
end

end

