function obj = trainSpBaTRIM(obj,theta_thilda, dX, varargin)

nn = 0;
for i =1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "min"
                nn = varargin{i+1};
        end
    end
end
trainRatio = obj.trainRatio;
obj.trainedalgorithm = "SpbTRIM";

% Training Data
[theta_thilda_Train,y_train,theta_thilda_test, y_test ] = splitData(theta_thilda, dX, trainRatio) ;
p = obj.maxpower;
obj.variablesVectors = string;

nmax = 20;
% Parameters
p = size(theta_thilda_Train, 2);
n_outputs = size(y_train, 2);
if size(theta_thilda, 2)>nmax

    knum = nmax;
else

    knum = size(theta_thilda, 2)-2;
end
% k-space
kvector       = round(linspace(nn, knum+1, knum));
loss_1        = zeros(knum-1 ,n_outputs );
spar          = zeros(knum-1 ,n_outputs );



lambda_0 = zeros(obj.mode.n, 1);
xi_0 = reduceXi( theta_thilda_Train,y_train,lambda_0);

% [theta_thilda_Train, theta_thilda_test]  = reduceTheta(theta_thilda_Train, theta_thilda_test, xi_0);

XI = containers.Map('KeyType', 'int32', 'ValueType', 'any');
for k = 1:length(kvector)


    % Optimizaiton objecitves 1

    for i=1:n_outputs
        [xi_trim, bayesParams]        = sparseBayesTrim3 (y_train(:, i) , theta_thilda_Train, kvector(k)) ;
        y_pred_1       = theta_thilda_test*xi_trim;
        y_true         = y_test(:, i);
        % Optimizaiton objecitves 2
        loss_1(k, i) =  sum((y_true - y_pred_1).^2, 'all');

        spar(k, :)   = nnz(xi_trim);

        mapVal(i).xi = xi_trim;
        mapVal(i).bayesParams = bayesParams;
    end

    XI(kvector(k)) = mapVal;

end

for i = 1:obj.mode.n
    f1         = loss_1(:, i);
    f2         = spar(:, i) ;%/max(spar(:, i));

    for j = 1:length(f1)
        RSS = f1(j);
        n_total = numel(y_true);
        k = f2(j);
        % Compute AIC and BIC
        AIC_vals(j) = 2*k + n_total*log(RSS / n_total);
        BIC_vals(j) = k*log(n_total) + n_total*log(RSS / n_total);

        gamma = 0.5;
        log_binom = gammaln(p+1) - gammaln(k+1) - gammaln(p - k + 1);
        % ebic = bic + 2 * gamma * log_binom;
        ebic(j) = BIC_vals(j)  + 2 * gamma * log_binom;
    end
    [~, maxIndex] = min(ebic);
    k_opt(i) = kvector(maxIndex);
end

for i = 1:obj.mode.n
    f1         = loss_1(:, i) /max(loss_1(:, i));
    f2         = spar(:, i) ;%/max(spar(:, i));
    if any(isnan(f1)) || any(isnan(f2))
        maxIndex = 1;
        f1  = loss_1(:, i);
        f2  = spar(:, i);
    elseif cov(f1)<1e-4

        maxIndex = 1;

    else

        maxIndex   =  traingleOPt(f1, f2 );
    end
    if obj.pltOption

        label = ["$\|Y -\Theta\Xi \|^2_2$";"$\|\xi_j\|_0$"];
        paretoPlot(f1, f2, maxIndex, ...
            "save", obj.savePlots, ...
            "var", obj.outputVars{i}, ...
            "noise", obj.noise, ...
            "label", label)

        plotObject(i).label = label;
        plotObject(i).save  = obj.savePlots;
        plotObject(i).var   = obj.outputVars{i};
        plotObject(i).f1 = f1;
        plotObject(i).f2 = f2;
        plotObject(i).kv = kvector;
    end

    k_opt(i) = kvector(maxIndex);
end
if obj.pltOption

    obj.plotObject = plotObject;

else
    obj.plotObject = []
end
p = size(theta_thilda_Train, 2);
x_opt = zeros(n_outputs, p);

for i=1:n_outputs
    mapOb = XI( k_opt(i));
    xi_trim = mapOb(i).xi;
    bayesParams = mapOb(i). bayesParams;


    spb{i} = bayesParams;
    x_opt(i, :) = xi_trim;
end

obj.Xi = x_opt';
obj.bayesParameters = spb;




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

function plotBayes(DIAGNOSTIC, HYPERPARAMETER,PARAMETER,noise, M, varargin)
label = ["\rho (density)"; "RICc"];
saveType = 0;
var      = "var";
plttitle = "SB";

filetitle = plttitle;

filetitle = strcat(filetitle,"-NL-", num2str(noise));
for i = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "save"
                saveType = varargin{i+1};
            case "var"
                var = varargin{i+1};
                filetitle = strcat(var,"-", plttitle);
                plttitle =  sprintf('\\bf{%s}\\rm\n%s', var, plttitle);
            case "noise"
                noiseLevel =  varargin{i+1};
                filetitle = strcat(plttitle,"-", num2str(noiseLevel));
        end
    end
end

w_infer						= zeros(M,1);
w_infer(PARAMETER.Relevant)	= PARAMETER.Value;

h1 = figure("Position",[100 100 800 400 ]);

subplot(2,2,1:2:3)
TITLE_SIZE	= 12;
lsteps	= length(DIAGNOSTIC.Likelihood);
plot(1:lsteps, DIAGNOSTIC.Likelihood, 'g-')
set(gca,'Xlim',[0 lsteps+1])
title('Log marginal likelihood trace','FontSize',TITLE_SIZE)

ax	= axis;
dx	= ax(2)-ax(1);
dy	= ax(4)-ax(3);
t_	= sprintf('Actual noise:   %.5f', noise);
text(ax(1)+0.1*dx,ax(3)+0.6*dy,t_,'FontName','Courier')
t_	= sprintf('Inferred noise: %.5f', 1/sqrt(HYPERPARAMETER.beta));
text(ax(1)+0.1*dx,ax(3)+0.5*dy,t_,'FontName','Courier')


% h2 = figure;

subplot(2,2,2)
h	= stem(w_infer,'filled');
set(h,'Markersize',3,'Color','r')
set(gca,'Xlim',[0 M+1])
t_	= sprintf('Inferred weights (%d)', length(PARAMETER.Relevant));
title(t_,'FontSize',TITLE_SIZE)
%
% Show the "well-determinedness" factors
%
subplot(2,2,4)
bar(DIAGNOSTIC.Gamma,'g')
axis([0 length(PARAMETER.Relevant)+1 0 1.1])
title('Well-determinedness (gamma)','FontSize',TITLE_SIZE)


if saveType
    parent_dir = strcat(pwd, '\Results\');
    folder_name = 'Pareto';

    full_path = fullfile(parent_dir, folder_name);

    if ~isfolder(full_path)

        mkdir(full_path)
    end

    filename_1 = strcat(filetitle, "log.png");
    filename_2 = strcat(filetitle, "weights.png");
    fullpath_1 = fullfile(full_path, filename_1);
    fullpath_2 = fullfile(full_path, filename_2);
    saveas(h1, fullpath_1)
    saveas(h2, fullpath_2)
end

end


