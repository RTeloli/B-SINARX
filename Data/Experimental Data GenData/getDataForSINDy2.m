function [X_train, y_train, testData, inputvariables, outputvariables, delays, arecur, vars] = getDataForSINDy2(nExperNeeded)


[z, u, v, x , a, t ] = getRandomized(nExperNeeded);%, "SignalType", "chirp");

figure
plot(x)
h = t(2)-t(1);
Y = [x  v  z  u  ]';
vars = {"x", "v", "z", "u"} ;
delays = [1 1 1 1] ;
arecur = logical([1 1 1  0 ]);
[input, output, ...
    inputvariables, outputvariables] = getTrainingData3(Y, ...
    "vars", vars, ...
    "recurIndex", arecur, ...
    "delays", delays, ...
    "noise", 0.000);


nExperNeeded           = 1;
excitationLevel        = "30";
Fileindex = 2;
SignalType             = "Chirp" ;
[ztest, utest, vtest, xtest , atest, t ] = getDataFolderFiles( ...
    "n_experiments",  nExperNeeded     ,...
    "amplitude",      excitationLevel,...
    "cuttoffFreq",    40,...
    "index",          Fileindex,...
    "signal",         SignalType );


Ytest = [xtest vtest ztest utest]';

[input_test, output_test]  = getTrainingData3(Ytest, ...
    "vars", vars, ...
    "recurIndex", arecur, ...
    "delays", delays);




X_train        = input' ;
y_train        = output';
%
testData.xtest = input_test';
testData.ytest = output_test';
end