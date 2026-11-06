function plotFreqFRF(y_pre, Ytest,h, leg, sveplt, frange) 
shift  = 0; 
  
 xtest   = Ytest(1, :);
utest   = Ytest(4, :); 
Fs = 1/h;
N      = length(xtest);
t_plot = 0:h:(N-1)*h; 
index  = t_plot > 0;  

t         = t_plot(index);
ytest     = Ytest(:, index);
Utest     = utest(index);
ytest     = ytest(:, 1:end-shift);
u1        = Utest(    1:end - shift);
t         = t(1:end-shift);

% Confidence interval again
Pc  = 99;
r_p = 0.5 * (100 + Pc);
r_m = 0.5 * (100 - Pc); 
x_p = cat(3, y_pre{:});  

 
clc
% Plot PSD for first variable
for i = 1:3
    xtrue    = ytest(i, :); 
    x_i_pred = squeeze(x_p(:, i, :))';

    x_i_pred = x_i_pred(:, index); 
    x_i_pred = x_i_pred(:, 1:end-shift);   
    
        frfPlotsDistribution(xtrue , x_i_pred,u1', Fs, ...
                         "freqrange", frange, ...
                         "legendInput", leg, ...
                         "snapsht", 1, ...
                         "ratio", 0.9, ...
                         "svplt", sveplt, ...
                         "isdb", 1, ...
                         "foldername", "Numerical");
    
end
end