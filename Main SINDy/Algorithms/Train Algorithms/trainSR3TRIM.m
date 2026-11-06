function obj = trainSR3TRIM(obj,theta_thilda, dX)


trainRatio = obj.trainRatio;
% Training Data
[theta_thilda_Train,y_train,theta_thilda_test, y_test ] = splitData(theta_thilda, dX, trainRatio) ;
p = obj.maxpower;
obj.variablesVectors = string;


% Parameters
p = size(theta_thilda_Train, 2);
n_outputs = size(y_train, 2);
if size(theta_thilda_Train, 2)>10

    knum = 10;
else

    knum = size(theta_thilda_Train, 2)-2;
end

% k-space
kvector       = round(linspace(2, knum+1, knum));
loss_1        = zeros(knum-1 ,n_outputs );
spar          = zeros(knum-1 ,n_outputs );



lambda_0 = zeros(obj.mode.n, 1);
xi_0 = reduceXi( theta_thilda_Train,y_train,lambda_0);

[theta_train_reduced, theta_test_reduced]  = reduceTheta(theta_thilda_Train, theta_thilda_test, xi_0);
p = size(theta_train_reduced, 2);

% Iteration for all k
j = 1;
 parfor k = 1:length(kvector)
    x_sol = zeros(n_outputs, p); 
    for i=1:n_outputs
        xi_trim        = trimmedSR3(theta_train_reduced, y_train(:, i), kvector(k));
        x_sol(i, :)    = xi_trim;
        % Test Data

        y_pred_1       = theta_test_reduced*xi_trim;
        y_true         = y_test(:, i);

        % Optimizaiton objecitves
        loss_1(k, i) = 0.5*norm(y_true - y_pred_1, 2);
        spar(k, i)   = calcDensity(xi_trim);
        Xi_col{i,k} = xi_trim';
    end

end

% [loss_1, kkk] = RICc( theta_thilda_test, y_test, Xi_col);

for i = 1:obj.mode.n
    f1         = loss_1(:, i)/max(loss_1(:, i));
    f2         = spar(:, i)/max(spar(:, i)); 
    if any(isnan(f1)) || any(isnan(f2))
            maxIndex = 1;
            f1  = loss_1(:, i);
            f2  = spar(:, i);
    else          
          
                maxIndex   =  traingleOPt(f1, f2 );
    end  
    if obj.pltOption

        label = ["$\|Y -\Theta\Xi \|^2_2$";"$\|\xi_j\|_0$"];
        paretoPlot(f1, f2, maxIndex, ...
            "title",  "TRIM", ...
            "save", obj.savePlots, ...
            "var", obj.outputVars{i}, ...
            "noise", obj.noise, ...
            "label", label)
    end

    k_opt(i) = kvector(maxIndex);
end

%%

p = size(theta_thilda_Train, 2);

x_opt = zeros(n_outputs, p);
x_opt2 = zeros(n_outputs, p);

for i=1:n_outputs
    xi_trim = trimmedSR3(theta_thilda_Train, y_train(:, i),k_opt(i));
    x_opt(i, :) = xi_trim;
    %
    % x_opt2(i, :) =  Xi_col{i, k_opt(i)};
end

obj.Xi = x_opt';




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



% 
% 
% 
% 
% 
% 
% 
% 
% i =2;
% clc; close all; 
% f1 =  (loss(:, i)) ;
% x1 =  (spar(:, i));
% x2 = lmb(:, i);
% re = isinf(x1);
% f1 = f1(~re);
% x1 = x1(~re);
% 
% gd = gradient(x1, f1) ; 
% [~, maxindex] = max(gradient(gd, x1));
% 
% gd2 = gradient(x1, f1); 
% [~, maxindexgd2] = max(gd2);
% 
% gd3 = gradient(x1, f1); 
% [~, maxindexgd3] = max(gd3); 
% 
% maxIndex2   =  traingleOPt(f1, x1 );
% 
% 
% figure
% subplot(2, 2, 1)
% plot(f1, x1)
% hold on
% plot(f1(maxindex) ,x1(maxindex)  , "or")
% plot(f1(maxIndex2),x1(maxIndex2), "ob")
% plot(f1(maxindexgd2),x1(maxindexgd2), "oblack")
% plot(f1(maxindexgd3),x1(maxindexgd3), "og")
% 
% 
% subplot(2, 2, 2)
% 
% plot(x1, gd)
% hold on
% plot(x1(maxindex),   gd(maxindex)  , "or")
% plot(x1(maxIndex2),  gd(maxIndex2), "ob")
% plot(x1(maxindexgd2),gd(maxindexgd2), "oblack")
% plot(x1(maxindexgd3),gd(maxindexgd3), "og")
% title("Gradient of loss against norm xi")
% 
% 
% subplot(2, 2, 3)
% plot(x1, gd2)
% hold on
% plot(x1(maxindex),    gd2(maxindex)  , "or")
% plot(x1(maxIndex2),   gd2(maxIndex2), "ob")
% plot(x1(maxindexgd2), gd2(maxindexgd2), "oblack")
% plot(x1(maxindexgd3), gd2(maxindexgd3), "og")
% 
% 
% subplot(2, 2, 4)
% plot(f1, gd3)
% hold on
% plot(f1(maxindex),    gd3(maxindex)  , "or")
% plot(f1(maxIndex2),   gd3(maxIndex2), "ob")
% plot(f1(maxindexgd2), gd3(maxindexgd2), "oblack")
% plot(f1(maxindexgd3), gd3(maxindexgd3), "og")
% title("Gradient of norm xi against loss ")
% 
% 
% 
% lmb_try = lmb(maxIndex2, :);
% xi_c = reduceXi(theta_train, dX_train,lmb_try);
% % xi_c(~(xi_c(:, i)==0),i)
% nnz(xi_c(:, i))
% 
% 
% lmb_try = lmb(maxindex, :);
% xi_c = reduceXi(theta_train, dX_train,lmb_try);
% % xi_c(~(xi_c(:, i)==0),i)
% nnz(xi_c(:, i))
% 
% lmb_try = lmb(maxindexgd2, :);
% xi_c = reduceXi(theta_train, dX_train,lmb_try);
% % xi_c(~(xi_c(:, i)==0),i)
% nnz(xi_c(:, i))
% 
% lmb_try = lmb(maxindexgd3, :);
% xi_c = reduceXi(theta_train, dX_train,lmb_try);
% % xi_c(~(xi_c(:, i)==0),i)
% nnz(xi_c(:, i))
% 
% 
% 
% 
% 
% 
% 
% function xi  = reduceXi(theta_x, y, lambda)
% [n, m] = size(theta_x);
% if n<m
%     theta_x=theta_x';
% end
% [n, m] = size(y);
% if n<m
%     y=y';
% end
% 
% n_basis = size(theta_x,2);
% % xi = (theta_x'*theta_x + 0.1*eye)\theta_x'*y;
% [~, n_outputs] = size(y) ;
% xi = zeros(n_basis, n_outputs);
% for i=1:n_outputs
%     [~, w] = sr3(theta_x,y(:, i),  "lam", lambda(i));
%     xi(:, i) = w; 
% end

% end