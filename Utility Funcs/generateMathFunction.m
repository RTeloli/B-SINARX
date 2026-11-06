function output = generateMathFunction(obj)
if any(strcmp(obj.stringExpNorm, "")) || any(isempty(obj.stringExpNorm))

    disp("Model not trianed")
else
    % h = obj.linearparams.h;
    % m = obj.linearparams.m;
    % c = obj.linearparams.c;
    % k = obj.linearparams.k;
    % output{1}.functions       =  @(x1,v1)h*v1+ x1;
    % output{1}.contained       =  "x1,v1";
    % output{1}.containedArray  =  logical([1 1 0 0]);
    % 
    % output{2}.functions       =  @(x1,v1,z1,u1) v1 +h/m*(u1 - (k*x1 + c*v1));
    % output{2}.contained       = "x1,v1,z1,u1";
    % output{2}.containedArray  = logical([1 1 1 1]);


    for i = 1:obj.mode.n
        [contained,stringExpression, containedArray] = getStringInvolvedMain(obj.vars, obj.stringExpNorm{i});
        % keyboard
        if isempty(contained)
            contained = "1";
            stringExpression = "1";  
        end
        inputStr  = char(sprintf('@(%s)', contained));
        stringExpression =  char(stringExpression);
        output{i}.functions      = str2func([inputStr, stringExpression]);
        output{i}.contained      = contained;
        output{i}.containedArray = logical(containedArray);


    end
end
end
