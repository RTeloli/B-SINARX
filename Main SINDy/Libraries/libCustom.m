function [theta, expression] = libCustom(x, expressions, vars)
n_expressions = length(expressions); 
[~, n] = size(x); 
theta = zeros(n_expressions, n); 

for i = 1:n_expressions
    exp   =  expressions{i}; 

    expression{i} = convertAbs(exp); 
    func = geneMathFuncs2(exp, vars); 
    X    = evaluateFunction(func, x); 
    theta(i, :) = X; 
end

end

function y = evaluateFunction(equations, A)

for i = 1:length(equations)
    a = A(equations{i}.containedArray,:); 
    b = num2cell(a,2);

    y(i, :) = equations{i}.functions(b{:});
end
end

function output = geneMathFuncs2(stringExpNorm, vars)
n = length(stringExpNorm);
if any(strcmp( stringExpNorm, "")) || any(isempty( stringExpNorm))

    disp("Model not trianed")
else 
    for i = 1:n
        [contained, stringExpression, containedArray] = generateCusFromExp( stringExpNorm{i}, vars);
       
        if isempty(contained)
            inputStr  = char(sprintf('@(%s)', "x"));
            stringExpression =  char("0*x");

            output{i}.functions      = str2func([inputStr, stringExpression]);
            output{i}.contained      = "x";
            output{i}.containedArray = logical(containedArray);
        else
            inputStr  = char(sprintf('@(%s)', contained));
            stringExpression =  strcat('real(',char(stringExpression), ')');
            output{i}.functions      = str2func([inputStr, stringExpression]);
            output{i}.contained      = contained;
            output{i}.containedArray = logical(containedArray);
        end

    end
end
end

function [containedVars, ex, contained] =  generateCusFromExp(expre, eqivars)
% vars = {'x', 'y'};
[vars, expre]  = makeVars(eqivars, expre);
operators = {'+', '-', '*', '/'};
genericFunctions = {"abs", "tan", "sin", "cos", "sinh", "cosh", "tanh", "exp", "max"};

clc
% expre     = 'x(y+x*8)';
expre     = strcat("1",expre);
expre     = char(expre);
n         = length(char(expre));
ex        = expre(1);

contained = cellfun(@(s) contains(expre, s), vars);

containedVars = strjoin([eqivars{contained}], ","); 
openBrace = false;
% ismember(expre, genericFunctions)
matches = cellfun(@(s) contains(expre, s), genericFunctions);
temp = [genericFunctions{matches}];
if isempty(temp)
    matched = '';
else
    matched = char(temp(end));
end
i = 1;
m = length(matched); 
while true
    % for i = 1:n-1
  
    term = expre(i); 
    nextTerm = expre(i+1); 
    
    if strcmp(expre(i), '(')
        openBrace = true;
    elseif strcmp(expre(i), ')')
        openBrace = false;
    end
    if sum(matches)>0 && m+i<n
        if strcmp(expre(i+1:m+i ), matched)
            newString = expre(i+1:m+i+1 );
            ex = strcat(ex, ".*",newString);

            i = i+m+1; 
            continue
        end
    end
    if strcmp(term, '(') || contains(term, operators )

        ex = strcat(ex,"__",nextTerm);
    elseif contains(nextTerm, operators )&&  ~isnan(str2double(term))

        ex = strcat(ex,"__",nextTerm);
    elseif contains(term, vars )    && ~contains(nextTerm, vars ) && openBrace
        
        ex = strcat(ex,nextTerm);
    elseif openBrace
        
        ex = strcat(ex,nextTerm);  
    elseif contains(term, '^') 
        
        ex = strcat(ex, "",nextTerm);
    elseif contains(nextTerm, '^')
         
        ex = strcat(ex, ".",nextTerm);
    else
        ex = strcat(ex, ".*","__", nextTerm);
    end



    i = 1 + i;
    if i >n-1
        break
    end
end 

for i = 1:length(eqivars)
    ex = strrep(ex, strcat('__', vars{i}),strcat(eqivars{i}));
    
end  

end

%%



