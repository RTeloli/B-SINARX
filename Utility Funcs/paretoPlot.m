function paretoPlot(f1, f2, maxIndex, varargin)
label = ["$\|Y -\Theta\Xi \|^2_2$";"$\|\xi_j\|_0$"];
saveType = 0;
var      = "var";
plttitle = "$L-Curve$";
filetitle = plttitle;
noiseLevel = 0;
for i = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "save"
                saveType = varargin{i+1};
            case "label"
                label = varargin{i+1};
            case "title"
                plttitle =  strcat(plttitle,"-", varargin{i+1});
            case "var"
                var = varargin{i+1};
                filetitle = strcat(var,"-", plttitle);
                plttitle =  sprintf('\\bf{%s}\\rm', var);
            case "noise"
                noiseLevel =  varargin{i+1};
        end
    end
end

filetitle = strcat(filetitle,"-", num2str(noiseLevel));

p3 = [f1(maxIndex) f2(maxIndex)];

figure
plot(f1, f2, "k","LineWidth",2)
hold on
grid on
plot(p3(1), p3(2), "o", "MarkerFaceColor", "red", "MarkerSize",8)
hold off
legend("Pareto Front", "Optimal", "interpreter", "latex")
% set(gca, "YScale", "log")
xlabel(label(1), "interpreter", "latex")
ylabel(label(2), "interpreter", "latex")
title(plttitle,  "interpreter", "latex")
set(gca, 'TickLabelInterpreter', 'latex')
set(gca,   'fontsize',18);
f1_rel = (max(f1)-min(f1))/max(f1)*100;
if f1_rel <2

    xlim([max(f1)*0.98 1.02*max(f1)])
end
if saveType

    parent_dir = strcat(pwd, '\Results\');
    folder_name = 'Pareto';

    full_path = fullfile(parent_dir, folder_name);

    if ~isfolder(full_path)

        mkdir(full_path)
    end

    filename = strcat(filetitle);
    fullpath = fullfile(full_path, filename);
    
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperSize', 0.5*[11 8.5]);
    % set(gcf, 'PaperPosition', [0 0 9.5 8.0]);
    set(gcf, 'PaperPosition', [0.25 0.25 10.5 8.0]); 
    print(gcf, '-dpdf', '-bestfit', fullpath,  '-vector');
end
end