%% =========== Numerical Case study - Bouc-Wen Hystersis ===================
% BSINARX Model Training and Testing Script
%  ------------------------------------------------------------------------
%  - Generates simulation data
%  - Prepares training and testing datasets
%  - Trains a sparse regression (SINDy) model with known physics
%  - Evaluates performance using Bayesian predictions
%  - Produces time response, hysteresis, and PSD plots
%  ========================================================================
clearvars; close all; clc
addpath(genpath('./'));

%% ========================================================================
%  Generate Data and Pre-Process
%  ========================================================================
nper = 0.0030;  % Noise level
[x, v, a, z, u, h] = getBoucWendata(nper, "tspan", [0 1]);    % Simulated data
Fs = 1/h;                                                     % Sampling frequency
Y  = [x; v; z; u']';                                          % Full dataset

% Variables for training
vars   = {"x", "v", "z", "u"};      % Input variables
delays = [1 1 1 1];                 % Number of delays per variable
arecur = logical([1 1 1 0]);        % Recursive inputs (1) vs exogenous (0)

% Prepare training data
[input, output, inputvars, outputvars] = getTrainingData3(Y, ...
    "vars", vars, ...
    "recurIndex", arecur, ...
    "delays", delays, ...
    "noise", nper);

% Known physics (linear system model)
[knownModel.xi_s, knownModel.expressions, knownModel.varibs] = BaselineDiscreteModel('numerical', 'linearbw',"timestep",h);

%% ========================================================================
%  Configure SINDy Model
%  ========================================================================
datamode.n     = sum(arecur);  % Number of recursive inputs
datamode.delay = delays;

% Library and training settings
libChoice   = "lib2Mod";   % Basis function library
power       = 2;           % Max polynomial power
trainRatio  = 0.80;        % Train/test split
n_train     = round(length(input) * trainRatio);

X_train     = input';
y_train     = output';
 

% Training algorithm
algoChoice = "SpbTRIM";   % (Check SndyModel.trainSpaseModel for options)

% Create SINDy model object
sndymodel = BSINARX(power,   datamode,libChoice, ...
    "inputvars", inputvars, ...
    "outvars", outputvars, ...
    "plotOpt", 1, ...
    "savePlots", 0, ...
    "trainType", algoChoice, ...
    "knownmodel", knownModel);

%% ========================================================================
%   Train Sparse Model
%  ========================================================================
disp("Training " + algoChoice + " ...")
sndyModelTrained02 = sndymodel.trainSpaseModel(X_train, y_train);
disp("Training Completed!")

%% ========================================================================
%   Display Trained Model Expressions
%  ========================================================================
clc
for i = 1:3
    fprintf("\n================= Normalized =================\n")
    fprintf('<strong>%s: f(X)</strong> (%s):\n\t%s\n', ...
        outputvars{i}, "SINA-2", sndyModelTrained02.stringExp{i})
end

%% ========================================================================
%  Save Trained Model (Optional, uncomment to save)
%  ========================================================================
% resp = input('<strong>Save trained model? (y/n)</strong>: ', 's');
% if contains(resp, "y")
%     FileName = algoChoice + "-" + num2str(round(nper*100)) + "-LR.mat";
%     folder   = "Trained Model";
%     save(fullfile(folder, FileName), "sndyModelTrained02");
% end

%% ========================================================================
%  Testing 
%  ========================================================================
type = "chirp";   % Test signal type ("chirp", "sin", "multisine")
[xtest, vtest, atest, ztest, utest, h, t_test, data] = ...
    getBoucWendata(0, "tspan", [0 1] ,"type", type);

Ytest = [xtest; vtest; ztest; utest'];

[input_test, output_test] = getTrainingData3(Ytest, ...
    "vars", vars, ...
    "recurIndex", arecur, ...
    "delays", delays);

%% ========================================================================
%  Bayesian Predictions
%  ========================================================================
models = sndyModelTrained02;
nsim   = 100;
ytest  = Ytest;
y_pre  = cell(nsim,1);
for i = 1:nsim
    [y_pre{i}, err] = models.bayesClosedLoop(ytest, delays, arecur);
end
ypre_nan = cellfun(@(s) any(any(isnan(s)')), y_pre);
y_pre(ypre_nan) =  []; 
%% ========================================================================
%  Time-Domain Plots with Confidence Intervals
%  ========================================================================
Pc = 99.9
sveplt = 0;     % Save plots (1/0)
plotLegend = {
    "True response"
    "BSINARX Corrected response"
    sprintf("$%.0f\\%%$ Confidence Interval", Pc)
    };

plotTime(y_pre, Ytest, h, plotLegend, sveplt, outputvars)

%% ========================================================================
%  Hystersis Plots
%  ========================================================================


sveplt = 0;
plotHysteris(y_pre, Ytest, plotLegend, sveplt)

%% ========================================================================
%  PSD Plots
% ========================================================================
sveplt = 0;
plotFreqFRF(y_pre, Ytest,h, plotLegend, sveplt)

