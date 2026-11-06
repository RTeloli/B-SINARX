function plotHysteris(y_pre, Ytest, leg, sveplt, ncycle, t_relPos) 
ztest   = Ytest(3, :);
utest   = Ytest(4, :); 
% ncycle = 1; 
cycle_indices = findCycle(utest', 1024, t_relPos, ncycle);

shift  = 1; 
range  = cycle_indices;    
zz     = ztest(1:end - shift);   
ytest  = Ytest(:, 1:end - shift);

% Confidence interval again
Pc  = 99;
r_p = 0.5 * (100 + Pc);
r_m = 0.5 * (100 - Pc); 
x_p = cat(3, y_pre{:}); 

% Predictions
x_i_pred = squeeze(x_p(1:end, 1, :))';
y1 = prctile(x_i_pred, r_p);
y2 = prctile(x_i_pred, r_m);
y3 = mean(x_i_pred, 1)';  

% Data for plotting
xtrue    = ytest(1, range);  
xdist.y1 = y1(range+shift); 
xdist.y2 = y2(range+shift); 
xdist.y3 = y3(range+shift)'; 
zzz      = zz(range); 

% Plot
plotterDistribution(zzz, xtrue, [], xdist, ...
    "legend"    , leg, ...
    "title"     , "Hysteresis plot", ...
    "savePlot"  , sveplt, ...
    "compare"   , "eearth2", ...
    "ratio"     , 0.1, ...
    "foldername", "Final", ...
    "type"      , "hystersis", ...
    "format"    , "pdf");  
end