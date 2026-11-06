function [y_pre, error]               = closedloopVarDelay1(obj, A, delays, arecur)
 
n = length(find(arecur));

[leng, m] = size(A);
if leng<m
    A=A';
    [leng, m] = size(A);
end 

maxdelay = max(delays); 

% Prediction function
f = @(x) predictFunction(obj, x);
 
% Variable Delay Mask
mask =  createDelayMask(delays, length(delays));

y_pre = zeros(leng , n); 

% Take the first step values from input upto max delay 
y_pre(1:maxdelay, :) = A(1: maxdelay, arecur);
ypast = y_pre( maxdelay, :);

% Error tracking for divergence
error  =   zeros(leng , n) ;
yTrue  = zeros(leng, n) ;
yTrue(1:2, :) = y_pre(1:2, :) ;
A0 =  A(1: maxdelay + 1-1, :)' ;
% Start loop from 
for i= 1:leng-maxdelay 
    clc
    % Input
    Adelay = A0 ;
    y_true = Adelay(arecur, end)';    % True value if correctly provided
     
    a0 = Adelay ;   
    a0 = a0.*mask; 

    a0(arecur,end ) = ypast';
    a0 = reshape(a0', 1, []);

    a0(isnan(a0)) =[] ;
     
    y = f(a0);

    y_pre(i+maxdelay, :) = y ;
    ypast = y;       
    A0 = circshift(A0,  [ 0 -1]); 
    A0(~arecur, :) =  A(i: maxdelay + i-1,~arecur)' ;

    error(i+maxdelay, :)  =  ( y_true- y );
    yTrue(i+maxdelay, :)  =   y_true ;
    
      clc
    checker = max(abs(y))/max(max(abs(A)));
    if any(any(isnan(error)))
         
         y_pre(i:end, :) =  NaN ; 

        fprintf("Relative error greater has Nan values! Divergent Solution!\n")
        return
    elseif checker> 20
        onesY = NaN*y_pre(i+maxdelay:end, :) ;
        y_pre(i+maxdelay:end, :) = onesY;
        fprintf("Relative error greater than 1! Divergent Solution!\n")
        return
    end 


end
end
function x = predFun2(x0, params)
xk   = x0(1);
vk_1 = x0(2);
vk   = x0(3); 
uk   = x0(4); 


m = params.m; 
c = params.c; 
k = params.k; 
h = params.h ;


ak = 1/m*(uk - c*vk - xk*k); 

x(1)  = xk + h*vk;
x(2)  = ak*h^2   + 2*vk - 1*vk_1; 
end

function mask = createDelayMask(delays, sz)
mask = ones( sz,max(delays));
[rowIdx, colIdx] = ndgrid(1:sz, 1:max(delays)  );
mask(  flip(colIdx > delays' , 2)) = NaN;


end