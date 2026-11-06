classdef BSINARX

    % This is class for B-SINARX object.
    % To initialize call (Naive):
    %            OBJ = BSINARX(LibraryFunctionName,MaxPower, DataMode)
    %    ** DataMode is an object/struct that must contain:
    %        DataMode.n     - Number of input 
    %        DataMode.delay - nu and ny, a vector of delays for each inputs [n_x1 n_x2, ..., n_u1]. 
    %                         It should be ordered the same way X        
    %  Class properties 
 

    properties
        library_function
        Xi
        stringExp
        stringExpNorm
        numTerms
        l
        maxpower
        mode
        variablesVectors = {}
        sparsity
        error
        funcs
        scaleMat
        vars
        outputVars
        bayesParameters
        noise
        pltOption
        savePlots
        trainedalgorithm
        MathFuncts
        modelKnowledge
        trainRatio
        plotObject
        CustomExp

    end

    methods
        function obj = BSINARX(power,mode,libFuncname ,  opts)
            arguments
                power        (1,1) double {mustBeInteger, mustBePositive}
                mode         struct
                libFuncname  (1,1) string {mustBeMember(libFuncname, ...
                    ["lib0","lib1","lib2","lib2Mod","ortho","sign","polysign"])} = "lib0"
                opts.OutVars         string = getStringCol(mode.n, "")
                opts.PlotOpt                = 0
                opts.Noise     (1,1) double {mustBeNonnegative} = 0
                opts.SavePlots (1,1) logical = false
                opts.InputVars       string = ""
                opts.TrainType (1,1) string {mustBeMember(opts.TrainType, ...
                    ["STLS","TRIM","SparseBayes","TSB","BSINARX","SR3","TRIMSR3"])} = "BSINARX"
                opts.KnownModel      struct = struct("xi_s", [], "expressions", [], "varibs", [])
            end

            obj.library_function = libFuncname;
            obj.maxpower         = power;
            obj.mode             = mode;
            obj.outputVars       = opts.OutVars;
            obj.pltOption        = opts.PlotOpt;
            obj.noise            = opts.Noise;
            obj.savePlots        = opts.SavePlots;
            obj.vars             = opts.InputVars;
            obj.trainedalgorithm = opts.TrainType;
            obj.modelKnowledge   = opts.KnownModel;

         
            if any(strlength(obj.vars) > 0)
                n = numel(obj.vars);
                x = 1000*rand(n, 100);
                    [~,  n_functions  ] = obj.libfunction(x , obj.maxpower); 
                disp(size(n_functions, 1))
            end
        end

        % Learning 
        function obj = trainSpaseModel(obj,x, dX)
            [n, m] = size(x);
            if n<m
                x=x';
            end
            [n, m] = size(dX);
            if n<m
                dX=dX';
            end

            [theta, obj.l, ~, functions ] = obj.libfunction(x', obj.maxpower);

            obj.funcs                          = functions;

            % Restructure
            xi_s        = obj.modelKnowledge.xi_s;
            expressions = obj.modelKnowledge.expressions;
            varibs      = obj.modelKnowledge.varibs;

            [theta_k, expressions]     = libCustom(x', expressions, varibs);
            y_know                     = theta_k'*xi_s;
            y_reformed                 = dX - y_know;


            % l2norm scale
            [theta_thilda, H] = l2scaleNorm(theta);

            obj.scaleMat   = H;
            theta_thilda   = theta_thilda';
            obj.trainRatio = 0.8;


            B = functions;
            S = [expressions{:}]'; 
            switch obj.trainedalgorithm
                % 1. Training STLS
                case "STLS"
                    obj = trainSTLS(obj,theta_thilda, y_reformed);
                case "TRIM"
                    % 2. Training TRIM optimization
                    [obj, toupdate] = trainTRIM(obj,theta_thilda, y_reformed);
                     % 
                     % obj.Xi(: ,~toupdate) = 0;
                case "SparseBayes"
                    % 3. Bayesian Sparse Regression Optimization
                    obj = trainSparseBayes(obj,theta_thilda, y_reformed);
                case "TSB"
                    % 4. Bayesian Thresholded Sparse Regression Optimization
                    obj = trainTSpB(obj,theta_thilda, y_reformed);
                case "BSINARX"
                    % 5. Training Trimmed Sparse Bayesian optimization
                    obj = trainSpBaTRIM(obj,theta_thilda, y_reformed, "min",2);
                case "SR3"
                    obj = trainSR3(obj,theta_thilda, y_reformed);
                case "TRIMSR3"
                    obj = trainSR3TRIM(obj,theta_thilda, y_reformed);
            end
 

            if ~isempty(obj.bayesParameters)
                for i = 1:n
                    % Bayesian Parameters
                    XiB = zeros(length(B), 1);
                    XiB(obj.bayesParameters{i}.index) = obj.bayesParameters{i}.param.Value;
                    xib = H\XiB;
                    xib = xib(obj.bayesParameters{i}.index);
                    xis = xi_s(:, i);
                    b = B(obj.bayesParameters{i}.index);
                    h = H(obj.bayesParameters{i}.index, obj.bayesParameters{i}.index);

                    Sigma1 = obj.bayesParameters{i}.Sigma;
                    [M , P, hUpdated]= createSimilarityMap(b, S, xis, xib, h);

                    a = M.values;
                    a = hUpdated*[a{:}]';
                    obj.bayesParameters{i}.param.Value = a;
                    sigma_m = P*Sigma1*P';

                    obj.bayesParameters{i}.Sigma       = sigma_m;   
                    terms =  M.keys;
                    obj.bayesParameters{i}.Expressions =  string(char([terms(:)]));
                    obj.bayesParameters{i}.H =  hUpdated;

                end
            else
                hUpdated = H;
            end

            m1_keys = B;
            m2_keys = S;

            M_keys = union(m1_keys, m2_keys, 'stable');

            P1 = build_embedding_matrix(m1_keys, M_keys);
            P2 = build_embedding_matrix(m2_keys, M_keys);
            for i = 1:n
                % Bayesian Parameters
                XiB = obj.Xi;
                xib = H\XiB(:, i);
                xis = xi_s(:, i);

                xi = P1*xib + P2*xis ;
                xi(abs(xi)<1e-6) = 0; 
                Xi_new(:, i) = xi;
                obj.scaleMat = hUpdated;

            end
            obj.Xi = Xi_new;
            % b = M.keys;
            obj.funcs    = M_keys;


            obj =  obj.makeEquationStructured();

        end        
      
        function [theta_X, l, string, funcs]  = libfunction(obj, x, p)
            funcs = "";
            switch obj.library_function
                case "lib0"
                    [theta_X, string,funcs, l]    = lib_0(x, p,obj.vars );
                case "lib1"
                    [theta_X, l, string,funcs]    = lib_1(x, p,obj.vars );

                    n = size(l, 2);
                    if isempty(obj.variablesVectors)
                        string = getStringCol(n, obj.vars);
                    end
                case "lib2"
                    n = size(x, 1);
                    col                    = getStringCol(n, obj.vars);
                    [theta_X, l ,string, funcs] = lib_2(x, p, col);

                case "lib2Mod"
                    n = size(x, 1);
                    col                    = getStringCol(n, obj.vars);
                    [theta_X, l ,string, funcs] = lib_2_Mod(x, p, col);
 
                case "ortho"
                    [theta_X, string, funcs, l] = libLegendre(x, p, obj.vars )  ;
                case "sign"
                         [theta_X, string, funcs, l] = libSign(x, p, obj.vars )  ;
             case "polysign" 
                         [theta_X, string, funcs, l] = libSignPoly(x, p, obj.vars )  ;
             
                otherwise
                    theta_X = x';
                    l = diag(ones(size(x, 1) ,1));
                    fprintf("Select valid library")

            end
             
 
        end

        function y_pre = openLoop(obj, A)
            theta_X = obj.libfunction(A', obj.maxpower);
            theta_X = obj.scaleMat\theta_X';

            theta_X = theta_X';c
            y_pre = theta_X*obj.Xi;
            y_pre = y_pre';
        end

        
        function [y_pre, eearth] = predictClosed(obj, A, varargin)
            [y_pre, eearth] = closedloop(obj,  A, obj.Xi);

        end

        function y_pre = predictClosed1(obj, A, delays,  arecur, varargin)
            y_pre = closedloopVarDelay1(obj, A, delays,  arecur);

        end
        function   [y_pre, errBay]  = bayesClosedLoop(obj, A, delays,  arecur, varargin)

            bayesobj = obj;
            n_ouput = sum(arecur);
            for i = 1:n_ouput
                mu     =  obj.bayesParameters{i}.param.Value ;
                Sigma  =   (abs(bayesobj.bayesParameters{i}.Sigma)).^2;

                xi                     = mvnrnd(mu', Sigma')';
                xi                     = obj.bayesParameters{i}.H\xi;
                terms                  = obj.bayesParameters{i}.Expressions;
                funcExperess           = GenerateFunction(terms, xi);
                bayesobj.MathFuncts(i) = geneMathFuncs2(funcExperess, obj.vars);
            end
            [y_pre, errBay]    = closedloopVarDelay3(bayesobj, A, delays,  arecur);

        end
        function  y_pre = predictClosed0(obj, A, delays, arecur)

            y_pre = closedloopVarDelay0(obj, A, delays, "",  arecur);

        end

        
        function obj =  makeEquationStructured(obj)
            if isempty(obj.funcs)|| strcmp(obj.funcs(end-1), "")
                if isempty(obj.variablesVectors)
                    n = size(obj.l, 2);
                    col = getStringCol(n, obj.vars);
                else
                    col = obj.variablesVectors;
                end
                obj.stringExp      = printFunctions(col, obj.l, obj.Xi);
                xi_norm            =  obj.Xi;
                obj.stringExpNorm  = printFunctions(col,  obj.l,xi_norm);
            else
                obj.stringExp      = GenerateFunction(obj.funcs,  obj.Xi);
                xi_norm            =  obj.Xi;
                obj.stringExpNorm  = GenerateFunction(obj.funcs,  xi_norm);
            end


            obj.MathFuncts = generateMathFunction(obj);

        end


    end
end
 
 
 
 