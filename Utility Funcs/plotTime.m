function plotTime(y_pre, Ytest, i_outout, h, leg, sveplt,outputvars) 

 xtest   = Ytest(1, :); 
snpsht = 1;     % Show snapshots (1/0)
shift  = 0;     
 
N      = length(xtest);
t_plot = 0:h:(N-1)*h;
index  = t_plot > 0;  

t      = t_plot(index);
ytest  = Ytest(:, index);
ytest  = ytest(:, 1:end-shift);
t      = t(1:end-shift);

% Confidence interval settings
Pc  = 99;                   % Confidence percentage
r_p = 0.5 * (100 + Pc); 
r_m = 0.5 * (100 - Pc); 
x_p = cat(3, y_pre{:});  

% Legend


% Unit scaling
unitMultiplier = [1000 1000 1];   % For nicer plots
 
    % Extract data
    xtrue = unitMultiplier(i_outout) * ytest(i_outout, :);
    x_i_pred = squeeze(x_p(:, i_outout, :))';
    % Confidence bounds
    y1 = unitMultiplier(i_outout) * prctile(x_i_pred, r_p);
    y2 = unitMultiplier(i_outout) * prctile(x_i_pred, r_m);
    y3 = unitMultiplier(i_outout) * mean(x_i_pred, 1)';

    y1 = y1(index) ;
    y2 = y2(index) ;
    y3 = y3(index)';


    % Shifted data 
    xdist.y1 = y1(1+shift:end) ;
    xdist.y2 = y2(1+shift:end) ;
    xdist.y3 = y3(1+shift:end);

    
    % Plot
    plotterDstrbtn(t, xtrue, [], xdist, ...
        "legend"    , leg, ...
        "title"     , outputvars{i_outout}, ...
        "savePlot"  , sveplt, ...
        "compare"   , "eearth2", ...
        "ratio"     , 0.95, ...
        "type"      , "time2", ...
        "foldername", "Final", ...
        "snapshot"  , snpsht, ...
        "format"    , "pdf");
 
end