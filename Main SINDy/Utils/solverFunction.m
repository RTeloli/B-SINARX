function y = solverFunction(obj, A)
h = obj.linearparams.h ; 
m = obj.linearparams.m ; 
c = obj.linearparams.c ; 
k = obj.linearparams.k ; 



equations  = obj.MathFuncts;
a = A(equations{3}.containedArray);
b = num2cell(a);
z = equations{3}.functions(b{:});

y(1) = A(1) + h*A(2); 
y(2) = A(2) + h*1/m*(A(4) - (k*A(1) + c*A(2) + A(3)));
y(3) = z; 



end

