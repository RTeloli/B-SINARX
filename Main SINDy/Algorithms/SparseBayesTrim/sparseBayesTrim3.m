function [Xi, bayesParams] = sparseBayesTrim3 (dXdt, Theta, k)
iterations = 1000;
OPTIONS		= SB2_UserOptions('iterations',iterations,...
    'diagnosticLevel', 1,...
    'monitor', 1000);

[M , N]= size(Theta');
n_outputs = size(dXdt, 2);
basisFunc = Theta;
Xi = zeros(M, n_outputs);

for ind = 1:n_outputs

    counter  = 1;
    try
        if k==0
                        bayesParams.hyper = 0;
                        bayesParams.param = 0;
                        bayesParams.diagn = 0;
                        bayesParams.index = 0;
                        bayesParams.Sigma = 0;
        else
            kk = k;
            % Xi =  xi0;
            target_train_i = dXdt(:, ind);
            basisFunc = Theta;
            sigma_2 = (std(target_train_i)^2*0.1)^0.5;

            orirgnalIndex = 1:M;

            SETTINGS	=  SB2_ParameterSettings('NoiseStd',sigma_2);
            % Initial estimate
            % lambda_0 = zeros(1, 1);
            % xi_0 = reduceXi( basisFunc,target_train_i,lambda_0);
            %
            % [param, HYPERPARAMETER, DIAGNOSTIC ]    = SparseBayes('Gaussian', basisFunc, target_train_i, OPTIONS, SETTINGS);
            % relev0    = param.Relevant;
            % relevVal0 = param.Value;
            %
            % ff = DIAGNOSTIC.Gamma;
            fmin =  0;
            ff = ones(length(orirgnalIndex), 1);
            relev0 = orirgnalIndex;
            relevVal0 = 0*orirgnalIndex;
            while true
                thresholder     =   ff>fmin;

                relev0   (~thresholder)    = [];

                basisFunc0       = Theta(:, relev0);
                if ~isempty(basisFunc0)
                    [param, HYPERPARAMETER, DIAGNOSTIC] = ...
                        SparseBayes('Gaussian', basisFunc0, target_train_i, OPTIONS, SETTINGS);

                else
                    param.Relevant = 0;
                    param.Value    = 0;

                end

                ff = DIAGNOSTIC.Gamma;
                fmin =  min(ff);

                relev0   = relev0(param.Relevant);
                if   counter >1000  ||  length(ff) <=kk
                    try
                        relevVal                = param.Value;
                        orirgnalIndex           = relev0;
                        Xi( relev0, ind) = relevVal;

                        sig2 = 1/HYPERPARAMETER.beta;
                        A = diag(HYPERPARAMETER.Alpha);
                        basisFunc       = Theta(:, orirgnalIndex);
                        Sigma = inv(HYPERPARAMETER.beta*(basisFunc'*basisFunc) +  diag(HYPERPARAMETER.Alpha));

                        bayesParams.hyper = HYPERPARAMETER;
                        bayesParams.param = param;
                        bayesParams.diagn = DIAGNOSTIC;
                        bayesParams.index = orirgnalIndex;
                        bayesParams.Sigma = Sigma;


                        break

                    catch e
                        keyboard
                    end
                else

                    counter = counter+1;

                    % continue
                end



            end
        end
    catch e
        disp(e.message)
        % keyboard
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
        % xi = (theta_x'*theta_x + 0.1*eye)\theta_x'*y;
        [~, n_outputs] = size(y) ;
        xi = zeros(n_basis, n_outputs);
        for i=1:n_outputs
            xi(:, i) = sparsifyDynamics(theta_x,y(:, i),  lambda(i), 1);
        end

    end


end