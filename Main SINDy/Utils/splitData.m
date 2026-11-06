function [xtrain,dX_train, xtest, dX_test, indices ] = splitData(x, y, trainRatio, varargin)

samplingFunc = "rand";
for i = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "samplfunc"
                samplingFunc = varargin{i+1};
        end
    end
end

switch samplingFunc
    case "direct"
        [n, ~] = size(y);
        trainLength = round(n*trainRatio);
        trainIndex  =  1:trainLength ;
        testInx     =  trainLength+1:n;
        
        indices.trainIndex = trainIndex  ;   
        indices.testInx    = testInx     ;   

        xtrain   = x(trainIndex, :);
        xtest    = x(testInx, :);
        dX_train = y(trainIndex, :);
        dX_test   = y(testInx, :);
    case "rand"
        [n, ~] = size(y);
        randIndx = randperm(n); 
        trainLength = round(n*trainRatio);
        trainIndex  =  randIndx(1:trainLength) ;
        testInx     =  randIndx(trainLength+1:n);

        indices.trainIndex = trainIndex  ;   
        indices.testInx    = testInx     ;   

        xtrain   = x(trainIndex, :);
        xtest    = x(testInx, :);
        dX_train = y(trainIndex, :);
        dX_test   = y(testInx, :);

end

