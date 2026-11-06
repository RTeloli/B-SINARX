function plotterDistribution(t, xtrue, xcompare, xdist, varargin)


%%%%%
%%% Plotter function
%%% Plots the variable vector provided in the columns of X matrix against t
%%% if options are provided defaults below will be chcanged
Pc = 99;

y1   = xdist.y1;
y2   = xdist.y2;
y3   = xdist.y3;
format = "png";
x = [xtrue;xcompare; y3];
[n, m] = size(x);
snapShot = 0;
savePlot = 0;
leg = getRandomLegend(n);
mode = "Response";
errcal = "eearth";
plotType = "time";
ratio = 0.1;
foldname = "";
figureC = 0;
errplt = 0; 
snapsveplt = 0;
for i = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}.lower()
            case "legend"
                leg = varargin{i+1} ;
            case "title"
                mode = varargin{i+1} ;
            case "snapsveplt"
                snapsveplt = varargin{i+1} ;

            case "saveplot"
                savePlot = varargin{i+1} ;
            case "snapshot"
                snapShot = varargin{i+1} ;
            case "compare"
                errcal = varargin{i+1} ;
            case "type"
                plotType= varargin{i+1} ;
            case "ratio"
                ratio = varargin{i+1} ;
            case "foldername"
                foldname = varargin{i+1} ;
            case "figure"
                figureC = varargin{i+1} ;
            case "interval"
                Pc = varargin{i+1};
            case "format"
                format = varargin{i+1};
            case "errplt"
                errplt = varargin{i+1};
        end
    end
end

switch plotType
    case "time"
        switch errcal
            case "eearth"
                accuracy = eearthCal(x);
                tit = sprintf(['\\bf{%s}\\rm\nEEARTH: ', ...
                    'Total = %.1f%%  ', ...
                    'Mag   = %.1f%%  ', ...
                    'Slope = %.1f%%  ', ...
                    'Phase = %.1f%%  '], ...
                    mode,...
                    accuracy.eet(1)*100, ...
                    accuracy.eem(1)*100, ...
                    accuracy.ees(1)*100, ...
                    accuracy.eep(1)*100);

            case "nise"
                accuracy = niseCal(x);
                tit = sprintf('\\bf{%s}\\rm\n', mode  );
                for i = 1:length(accuracy)
                    tit = [tit sprintf('NISE_{%s} = %.2f%%  ', leg{i+1}, accuracy(i)*100 )];
                end
            case "eearth2"

                accuracy = eearthTotal(x);
                tit = sprintf('\\bf{%s}\\rm\n', mode  );
                for i = 1:length(accuracy)
                    tit = [tit sprintf('EEA_{%s} = %.2f%%  ', leg{i+1}, accuracy(i)*100 )];
                end

            otherwise
                tit = "";
        end
        for i = 2:n
            error_displacement(i-1, :) = x(1, :) - x(i, :);
        end

        [n, m] = size(x);
        color = [];
        % color = rand(3,n);
        % color(:, 1) = 1;
        % color(3,:) = 0;
        color = [color, "black", "red", "blue", "green", "cyan", "magenta"];
        lWidth = [ 5 1.5 1.5 1.5, 1.5, 1.5];
        % lWidth = [ 5 1 1 1, 1, 1 ];
        lineType = [{'.'}, {'-'}, {'-'}, {'-'}, {'-'}, {'-'}];


        if figureC==0
            figureC = figure;
            set(figureC, "Position",[50 50 1000 500])
        end

        subplot(4, 1, 1:3)

        for i = 1:n
            plot(t, x(i, :) , lineType{i}, "Color" ,color(:,i), "LineWidth",lWidth(i))
            hold on
        end
        patch([t fliplr(t)]'  , [y1 fliplr(y2)] , 'g',  'FaceAlpha', 0.5, 'EdgeColor', 'none');


        xlim([min(t) max(t)])
        switch mode
            case "x"
                unit = "(m)";
            case "v"
                unit = "(m/s)";
            case "z"
                unit = "(N)";
            otherwise
                unit = "[unit]";
        end
        ylab = sprintf('${%s}$', mode);
        ylabel(strcat(ylab, unit),'Interpreter', 'latex', 'FontWeight', 'bold');


        grid on
        hold off
        sgtitle(tit, 'Interpreter', 'tex')
        xlabel("t(s)")
        xlim([min(t) max(t)])
        legend(leg)

        subplot(4, 1, 4)
        plot(t, error_displacement, "LineWidth",1.5)
        legend(leg(2:end))


        xlim([min(t) max(t)])
        xlabel("t(s)")
        grid on

        ylabel("Error")

    case "hystersis"
        zzz = t;
        y_plot_true = xtrue;


        xx1   = xdist.y1;
        xx0   = xdist.y2;
        xx4   = xdist.y3;
        if ~isempty(xcompare)
            x_plot = [xx4; xcompare];
        else
            x_plot =  xx4 ;
        end
        color  = [];
        color  = [color, "red", "blue", "green", "cyan", "magenta"];

        lWidth = [   1.5 1.5 1.5, 1.5, 1.5];

        figure
        plot(y_plot_true, zzz , "oblack", "MarkerFaceColor","black", "LineWidth",lWidth(1))
        hold on
        for i=1:size(x_plot, 1)
            plot(x_plot(i, :), zzz , "Color",color(i),"LineWidth",lWidth(i+1))
            hold on
        end
        zzz = zzz(:);      % vertical axis (Y)
        xx0 = xx0(:);      % lower X bound
        xx1 = xx1(:);      % upper X bound

        % Construct closed polygon
        x_fill = [ xx0; flipud(xx1); ];  % left side up, right side down
        y_fill = [ zzz; flipud(zzz); ];  % zzz up, zzz down

        % Plot patch
        fill(x_fill, y_fill, 'k',  'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold off
        legend(leg,  "location", "best",  "Interpreter","latex")

        xlabel("Displacement(m)", "Interpreter","latex")
        ylabel("Hystersis Force (N)", "Interpreter","latex")
        grid on

        set(gcf,'units','inches','renderer', 'painters');
        set(gca,'TickLabelInterpreter','Latex','fontsize',16);

    case "time2"
        switch errcal
            case "eearth"
                accuracy = eearthCal(x);
                tit = sprintf(['\\bf{%s}\\rm\nEEARTH: ', ...
                    'Total = %.1f%%  ', ...
                    'Mag   = %.1f%%  ', ...
                    'Slope = %.1f%%  ', ...
                    'Phase = %.1f%%  '], ...
                    mode,...
                    accuracy.eet(1)*100, ...
                    accuracy.eem(1)*100, ...
                    accuracy.ees(1)*100, ...
                    accuracy.eep(1)*100);

            case "nise"
                accuracy = niseCal(x);
                tit = sprintf('\\bf{%s}\\rm\n', mode  );
                for i = 1:length(accuracy)
                    tit = [tit sprintf('NISE_{%s} = %.2f%%  ', leg{i+1}, accuracy(i)*100 )];
                end
            case "eearth2"

                % accuracy = eearthTotal(x);
                tit = sprintf('\\bf{%s}\\rm\n', mode  );
                % for i = 1:length(accuracy)
                %     tit = [tit sprintf('EEA_{%s} = %.2f%%  ', leg{i+1}, accuracy(i)*100 )];
                % end

            otherwise
                tit = "";
        end
        for i = 2:n
            error_displacement(i-1, :) = x(1, :) - x(i, :);
        end

        [n, m] = size(x);
        color = [];
        color = [color, "black", "red", "blue", "green", "cyan", "magenta"];
        lWidth = [ 5 1.5 1.5 1.5, 1.5, 1.5];
        lineType = [{'.'}, {'-'}, {'-'}, {'-'}, {'-'}, {'-'}];


        if figureC == 0
            figureC = figure;
            set(figureC, "Position",[50 50 1000 500])
        end


        for i = 1:n
            plot(t, x(i, :) , lineType{i}, "Color" ,color(:,i), "LineWidth",lWidth(i))
            hold on
        end
        patch([t fliplr(t)]'  , [y1 fliplr(y2)] , 'k',  'FaceAlpha', 0.2, 'EdgeColor', 'none');


        xlim([min(t) max(t)])
        switch mode
            case "x"
                unit = "[mm]";
            case "v"
                unit = "[mm/s]";
            case "z"
                unit = "[N/kg]";
            otherwise
                unit = "[unit]";
        end
        ylab = sprintf('${%s}$', mode);
        ylabel(strcat(ylab, unit),'Interpreter', 'latex', 'FontWeight', 'bold');


        grid on
        sgtitle(tit, 'Interpreter', 'latex')
        xlabel("t(s)")
        xlim([min(t) max(t)])

         if errplt
        plot(t, error_displacement, "LineWidth",1.5)

        leg{end+1} = "Error";
        
         end
         legend(leg, 'Interpreter', 'latex')
         xlim([min(t) max(t)])
         ylim([1.6*min(y3) 2*max(y3)])
         xlabel("t(s)",'Interpreter', 'latex', 'FontWeight', 'bold');
         grid on
         set(gcf,'units','inches','renderer', 'painters');
         set(gca,'TickLabelInterpreter','Latex','fontsize',18);
        % exportgraphics(gca,['Loss_Evolution_Er_Variable.eps'],'BackgroundColor','none')

end
if savePlot
    parent_dir = strcat(pwd, '\Results\',foldname);
    folder_name = plotType;

    full_path = fullfile(parent_dir, folder_name);

    if ~isfolder(full_path)

        mkdir(full_path)
    end
    unique_id = string(datetime('now','TimeZone','local','Format','HHmmss'));
    
    if snapsveplt
        fileName = strcat( "Plot-", mode, "-", num2str(unique_id), "SnapShot.", format);
    else
        fileName = strcat( "Plot-", mode, "-", num2str(unique_id), ".", format);
    end
    

    folderpath = fullfile(full_path, fileName);
    switch format
        case "pdf"
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperOrientation', 'landscape');
            set(gcf, 'PaperSize',0.7*[11 8.5]); 
            print(gcf, '-dpdf', '-fillpage', folderpath,  '-vector');
           otherwise
            saveas(gcf, folderpath)
    end
end
if snapShot
    range         = round(ratio*length(t));
    randRange     = randperm(length(t)-range);
    mode_snap     = strcat(mode);
    range         = randRange:(randRange+range);
    xdistRange.y1 = xdist.y1(range);
    xdistRange.y2 = xdist.y2(range);
    xdistRange.y3 = xdist.y3(range);
    if ~isempty(xcompare)
        xcom = xcompare( range);
    else
        xcom = [];
    end
    plotterDistribution(t(range), xtrue(:, range),xcom,xdistRange, ...
        "legend"  ,   leg , ...
        "title"   ,    mode_snap  , ...
        "savePlot", savePlot, ...
        "snapshot", 0, ...
        "compare"  , errcal, ...
        "foldername", foldname, ...
        "type", plotType, ...
        "format", format, ...
        "snapsveplt", 1);
end
end

function eearth = eearthCal(x)
[n, m] = size(x);
for i = 1:n-1
    if any(isnan(x(i, :))) || any(isnan(x(i+1, :)))
        eearth.eet(i) = 0;
        eearth.eep(i) = 0;
        eearth.eem(i) = 0;
        eearth.ees(i) = 0;

    else
        [   eearth.eet(i), ...
            eearth.eep(i), ...
            eearth.eem(i), ...
            eearth.ees(i)] = EEARTH2(x(i, :)', x(i+1, :)');
    end
end
end
function eearth = eearthTotal(x)
[n, m] = size(x);
for i = 1:n-1
    if  any(isnan(x(i+1, :)))
        eearth(i) = 0;


    else
        eearth(i)  = EEARTH2(x(1, :)', x(i+1, :)');
    end
end
end

function acc = niseCal(x)
[n, m] = size(x);
for i = 1:n-1
    if any(isnan(x(i+1, :)))
        acc(i) = 0;

    else
        acc(i) = nise(x(1, :)', x(i+1, :)');
    end
end
    function error                  = nise(a, b)
        phi_aa = sum(a.^2)/length(a);
        phi_bb = sum(b.^2)/length(b);
        phi_ab = sum(b.*a)/length(b);

        error = 2*phi_ab/(phi_aa+phi_bb);
    end


end



function leg = getRandomLegend(n)
defaultText = "Response";
for i = 1:n
    leg{i} = strcat(defaultText, "_", num2str(i));
end
end




