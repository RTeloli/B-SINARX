function [input, output, inputvariables, outputvariables] =  getTrainingData3(Y, varargin)

[n, m]  = size(Y);
if n>m
    Y = Y'; 
    [n, m] = size(Y); 
end

vars = getStringCol(n,"");
recurVars = logical([0 0 1 0  ]); 
delays = ones(n, 1); 
 nper = 0.00; 
for i = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "vars"
                vars = varargin{i+1}; 
            case "recurIndex"
                recurVars = varargin{i+1}; 
            case "delays"
                delays = varargin{i+1};
            case "noise"
                nper = varargin{i+1};
        end
    end
end
for i = 1:n
    Y(i, :) = addnoise(Y(i, :), nper)'; 
end
maxdelay = max(delays); 

 

 
input = [];

for i = 1:length(delays)
    for j = delays(i):-1:1
        input  = [input; Y(i, j:end-maxdelay+j-1)];
    end
end 
output = Y(recurVars, max(delays(recurVars))+1:end)';  
input = input'; 
  
inputvariables  =  getInputString(vars, delays);
outputvariables = vars(recurVars); 
             
end

 
function outputStng= getInputString(inputSting, delay)
outputStng = {};
n_input = length(inputSting);
l = 1;
for i = 1:n_input
    for j =   delay(i):-1:1
        outputStng{l} = strcat(inputSting{i}, num2str(j));
        l = l+1;
    end
end
end