function obj = trainTSpB(obj ,theta_norm, dX)
obj.trainedalgorithm = "TSB";

trainRatio = obj.trainRatio;

lambdan = vecnorm(dX, 2);
if   any(abs(lambdan)<eps)
    lambdan = eps+lambdan;
end


[theta_train,dX_train,theta_test, dX_test ] = splitData(theta_norm, dX, trainRatio) ;

lambda_0            = zeros(obj.mode.n, 1);
xi_0                = bayesianSINDy(dX_train, theta_train,lambda_0);

[theta_train_reduced, theta_test_reduced]  = reduceTheta(theta_train, theta_test, xi_0);

% xi_0_cpy            = xi_0;
% [~,inx]             = maxk(xi_0, 3);

% for i = 1:obj.mode.n
    % xi_0_cpy(inx(:, i), i) =NaN;
% end
% xi_0_cpy(xi_0_cpy ==0) = NaN;
n = 5000;
% lmbmin = 1*min(abs(xi_0_cpy), [], 1, 'omitnan');
% lmbmax = 1*max(abs(xi_0_cpy), [], 1, 'omitnan');
% 
% lmb =   lmbmin' + (lmbmax -lmbmin )'*linspace(0, 1, n );
% lmb = lmb';

lambda_rel = logspace(-3, 1, n );
lmb = lambda_rel'*lambdan;

spar = zeros(n, size(lmb, 2));

nnormxi = zeros(n, size(lmb, 2));
clc
fprintf("Running Optimization\n")

parfor i = 1:n
    xi{i}      = bayesianSINDy( dX_train,theta_train_reduced, lmb(i, :));
    spar(i, :) = size(xi{i}, 1) - sum(xi{i}==0);
    nnormxi(i, :) = vecnorm(xi{i}, 1); 
end
[~, kk, loss] = RICc(theta_test_reduced, dX_test , xi);


for i = 1:obj.mode.n

    above2 = (spar(:, i)>1);
    f1        =  loss(above2, i)/max(loss(above2, i));
    f2        =  spar(above2, i)/max(spar(above2, i));
    if any(isnan(f1)) || any(isnan(f2))
            maxIndex = 1;
            f1  = loss(:, i);
            f2  = spar(:, i);
    else          
          
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

xi_c = bayesianSINDy(dX_train,theta_train,lmb_try);

obj.Xi = xi_c;

end
