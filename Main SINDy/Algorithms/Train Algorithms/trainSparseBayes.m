function obj = trainSparseBayes(obj,theta_thilda, dX)

obj.trainedalgorithm = "SparseBayes";

basis_train  = theta_thilda;
n_outputs    =  obj.mode.n ;

% Training Data
xi = zeros(n_outputs, size(basis_train, 2));
xiimportvar = zeros(n_outputs, size(basis_train, 2));
% xi_0 = basis_train\dX;
clc

for i = 1:n_outputs
    target_train_i = dX(:, i);
    % SBS from demo
    iterations	= 5000;
    OPTIONS		= SB2_UserOptions('iterations',iterations,...
        'diagnosticLevel', 1,...
        'monitor', 1000);
    SETTINGS	= SB2_ParameterSettings('NoiseStd',1000);

    [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = ...
        SparseBayes('Gaussian', basis_train, target_train_i, OPTIONS, SETTINGS);
    xi(i, PARAMETER.Relevant) = PARAMETER.Value;

    if obj.pltOption
        plotBayes(DIAGNOSTIC, ...
            HYPERPARAMETER, ...
            PARAMETER, ...
            0, ...
            size(basis_train, 2), ...
            "save", obj.savePlots, ...
            "var", obj.outputVars{i})
    end

    var = 1./HYPERPARAMETER.Alpha;
    xiimportvar(i, PARAMETER.Relevant)  = var ;
    basisFunc = basis_train(:, PARAMETER.Relevant);
    Sigma = inv(HYPERPARAMETER.beta*(basisFunc'*basisFunc) +  diag(HYPERPARAMETER.Alpha));

    obj.bayesParameters{i}.Sigma  = Sigma;
    obj.bayesParameters{i}.hyper  = HYPERPARAMETER;
    obj.bayesParameters{i}.param  = PARAMETER;
    obj.bayesParameters{i}.diagn  = DIAGNOSTIC;
    obj.bayesParameters{i}.index  = PARAMETER.Relevant;


    obj.Xi = xi';

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

h1 = figure;
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


h2 = figure;

subplot(2,1,1)
h	= stem(w_infer,'filled');
set(h,'Markersize',3,'Color','r')
set(gca,'Xlim',[0 M+1])
t_	= sprintf('Inferred weights (%d)', length(PARAMETER.Relevant));
title(t_,'FontSize',TITLE_SIZE)
%
% Show the "well-determinedness" factors
%
subplot(2,1,2)
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


