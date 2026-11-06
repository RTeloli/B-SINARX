function Xi= bayesianSINDy(dXdt, Theta, lambda)
iterations = 1000;
OPTIONS		= SB2_UserOptions('iterations',iterations,...
    'diagnosticLevel', 1,...
    'monitor', 1000);

[M , N]= size(Theta');
n_outputs = size(dXdt, 2);
basisFunc = Theta;
k  = 1;
Xi = zeros(M, n_outputs);
for ind = 1:n_outputs 
 
    % Xi =  xi0; 
    target_train_i = dXdt(:, ind);
    basisFunc = Theta;
    sigma_2 = (std(target_train_i)^2*0.1)^0.5;
    SETTINGS	=  SB2_ParameterSettings('NoiseStd',sigma_2);
    % Initial estimate
    param =  SparseBayes('Gaussian', basisFunc, target_train_i, OPTIONS, SETTINGS);
    relev0 = param.Relevant;
    relevVal0 = param.Value;
     
    orirgnalIndex = 1:M; 
   while true
        thresholder     =  ~(abs(relevVal0)>lambda(ind));    
         
        % Xi(thresholder) = [];               
        relev0(thresholder) = [];
        relevVal0(thresholder) = [];
  
        basisFunc       = basisFunc(:, relev0);
        orirgnalIndex   = orirgnalIndex(relev0);

        if  ~any(thresholder) || isempty(basisFunc) |  k >100
 
            Xi( orirgnalIndex, ind) = relevVal0;
            break
        end
        [param, HYPERPARAMETER, DIAGNOSTIC] = ...
            SparseBayes('Gaussian', basisFunc, target_train_i, OPTIONS, SETTINGS);
         k = k+1;


        relev0 = param.Relevant;
        relevVal0 = param.Value;

       
    end


    
end