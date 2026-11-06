function y = predictFunction(obj, A)
equations = obj.MathFuncts; 

for i = 1:length(equations)
    a = A(equations{i}.containedArray);
    b = num2cell(a);
    y(i) = equations{i}.functions(b{:});
end
end
