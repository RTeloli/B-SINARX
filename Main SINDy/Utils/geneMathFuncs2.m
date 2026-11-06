function output = geneMathFuncs2(stringExpNorm, vars)
n = length(stringExpNorm);
if any(strcmp( stringExpNorm, "")) || any(isempty( stringExpNorm))

    disp("Model not trianed")
else 
    for i = 1:n
        [contained,stringExpression, containedArray] = getStringInvolvedMain( vars,  stringExpNorm{i});
       
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
