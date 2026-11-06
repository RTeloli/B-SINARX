clearvars; close all; clc
addpath(genpath('./'));
%% Generate Data
nExperNeeded         = 2 ;
[z, u, v, x , a, t ] = getRandomized(nExperNeeded);

nper = 0.00;
h = t(2)-t(1);
Y = [x  v  z  u  ]';
vars = {"x", "v", "z", "u"} ;
delays = [1 1 1 1]   ;
arecur = logical([1 1 1  0 ]);
[input, output, ...
    inputvars, outputvars] = getTrainingData3(Y, ...
    "vars",       vars, ...
    "recurIndex", arecur, ...
    "delays",     delays, ...
    "noise",      nper);

figure
plot(x)
hold on
plot(input(:, 1), "--")
%%
datamode.n     = sum(arecur);  % Number of recursive inputs
datamode.delay = delays;
datamode.h     = h;            % Time step


% Library and training settings
libChoice   = "lib0";   % Basis function library
power       = 3;           % Max polynomial power
trainRatio  = 0.90;        % Train/test split
n_train     = round(length(input) * trainRatio);

X_train = input';
y_train = output';

[xi_s, expressions, varibs] = BaselineDiscreteModel('experimental','linear_ex', "timestep", h);

% Known physics model
knownModel.xi_s        = xi_s;
knownModel.expressions = expressions;
knownModel.varibs      = varibs;

% Training algorithm
algoChoice = "SpbTRIM";   % (Check BSINARX.trainSpaseModel() for options)

 

%  SINDy Model

% Create SINDy model object
bsinarxmodel =  BSINARX(power,   datamode,libChoice, ...
    "InputVars", inputvars, ...
    "OutVars",   outputvars, ...  
    "TrainType", "BSINARX", ...
    "knownmodel", knownModel);

Fs = 1/h;

 
%% ========================= Train Sparse Model ===========================  
%  ========================================================================

disp(strcat("Training ", algoChoice))
bsinarxmodel    = bsinarxmodel.trainSpaseModel(X_train, y_train);


%% =============== Display Trained Model Expressions=======================
%  ========================================================================
clc
for i = 1:datamode.n
    fprintf("\n================= Normalized =================\n")
    fprintf('<strong>%s: f(X)</strong> (%s):\n\t%s\n', ...
        outputvars{i}, "SINA-2", bsinarxmodel.stringExp{i})
end
 
 
%%
nExperNeeded           = 1;
excitationLevel        = "15";
Fileindex              = 2;
SignalType             = "Chirp" ;
[ztest, utest, vtest, xtest , atest, t ] = getDataFolderFiles( ...
    "n_experiments",  nExperNeeded     ,...
    "amplitude",      excitationLevel,...
    "cuttoffFreq",    80,...
    "index",          Fileindex,...
    "signal",         SignalType );

Ytest = [xtest vtest ztest utest]';

[input_test, output_test]  = getTrainingData3(Ytest, ...
    "vars", vars, ...
    "recurIndex", arecur, ...
    "delays", delays);

%%
nsim  = 50;
ytest = Ytest ;
y_pre = cell(nsim,1);
for i = 1:nsim
    y_pre{i}  = bsinarxmodel.bayesClosedLoop(ytest, delays,  arecur) ;
end

ypre_nan = cellfun(@(s) any(any(isnan(s)')), y_pre);
y_pre(ypre_nan) =  []; 
%% ========================================================================
%  Time-Domain Plots with Confidence Intervals
%  ========================================================================

Pc = 95;
sveplt = 0;     % Save plots (1/0)
plotLegend = {
    "True response"
    "BSINARX Corrected response"
    sprintf("$%.0f\\%%$ Confidence Interval", Pc)
    };
for i_outout = 1:datamode.n
    plotTime(y_pre , Ytest, i_outout, h, plotLegend, sveplt, outputvars)
end
 
%% ========================================================================
%  Hystersis Plots
%  ========================================================================

t_relPos = 0.22;
ncycle   = 1.5;
sveplt   = 0;

plotHysteris(y_pre, Ytest, leg, sveplt, ncycle, t_relPos)


 




