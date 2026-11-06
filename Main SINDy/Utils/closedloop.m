function [y_pre, error]               = closedloop(obj,  A, xi)
n = obj.mode.n;
[leng, m] = size(A);
if leng<m
    A=A';
    [leng, m] = size(A);
end
delay = obj.mode.delay; 
f = @(x) predictFunction(obj, x);
 

error = [];
y_pre = [];
yTrue = [];

a0 = A(1,:);
y = f(a0);

error  = [error; diff([y; y])];
yTrue = [yTrue y'];
y_pre = [y_pre y'];
for i= 2:leng
    clc
    fprintf("<strong>%s\nClosed loop prediction</strong>: %.2f%%\n", obj.trainedalgorithm ,i*100/leng)
    x_0 = reshape(a0, [delay m/delay])';
    x_0 = circshift(x_0, [0 -1]);

    x_unreplaced      = reshape(A(i,:), [delay m/delay])';
    y_true            = x_unreplaced(1:n, end);
    x_0(n+1:end, end) = x_unreplaced(n+1:end, end);
    x_0(1:n, end)     = y;

    a0 = reshape(x_0', 1, []);

    y = f(a0);

    error  = [error; diff([y_true'; y])];
    yTrue  = [yTrue y_true];
    
clc
    checker = max(abs(y))/max(max(abs(A)));
    if any(any(isnan(error)))
        onesY = NaN*y_pre(:, end).*ones(n, size(A, 1)) ;
        y_pre = [y_pre onesY(:, i+1:end)];

        fprintf("Relative error greater has Nan values! Divergent Solution!\n")
        return
    elseif checker> 2
        onesY = NaN*y_pre(:, end).*ones(n, size(A, 1)) ;
        y_pre = [y_pre onesY(:, i:end)];
        fprintf("Relative error greater than 1! Divergent Solution!\n")
        return
    end
    y_pre = [y_pre y'];


end
end


