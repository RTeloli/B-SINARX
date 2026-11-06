function frfPlotsDistribution(xtrue, xdist,xexcit, fs, varargin)
% Plot PSD of signals with distribution (mean + confidence interval) in dB
% xtrue: true signals (nsignal x L)
% xdist: distribution of signals (ndist x L)
% fs: sampling frequency
% Optional inputs:
%   'freqrange', [fmin fmax]
%   'legendinput', {'label1','label2',...}

% Parse varargin
frange = [];
leg    = [];
i      =  1;
snapShot = 0;
savePlot = 0;
N = length(xexcit);
ratio = 0.95;
foldname = "FRF";
isdb = 0;
format = "pdf";

foldnameSpecial = ""; 
% PSD parameters 
Pc       = 99;                  % confidence interval
r_p      = 0.5*(100 + Pc);
r_m      = 0.5*(100 - Pc);
scale    = "linear";
nfs      = fs;
wind = kaiser(N,1);


while i <= length(varargin)
    if ischar(varargin{i}) || isstring(varargin{i})
        switch lower(string(varargin{i}))
            case "freqrange"
                frange = varargin{i+1};
                i = i + 1;
            case "legendinput"
                leg = varargin{i+1};
                i = i + 1;
            case "snapsht"
                snapShot = varargin{i+1};
            case "ratio"
                ratio = varargin{i+1};
            case "svplt"
                savePlot= varargin{i+1};
            case "isdb"
                isdb= varargin{i+1};
            case "foldername"
                foldnameSpecial = varargin{i+1};
        end
    end
    i = i + 1;
end



% Compute PSDs for distribution
[nDist, ~] = size(xdist);
Hdsit    = {};

for i = 1:nDist
    [Hxy,f] = tfestimate(xexcit,xdist(i,:), wind, [], nfs, fs);
    if isdb
        Hdsit{i} = mag2db(abs(Hxy));
    else
        Hdsit{i} =  abs(Hxy);

        scale = "log";
    end
end

% Compute mean and confidence intervals in dB
matHxy =  cell2mat(Hdsit);
yMean_dB = mean(matHxy, 2);% 10*log10
yHigh_dB = (prctile(matHxy, r_p, 2));% 10*log10
yLow_dB  = (prctile(matHxy, r_m, 2));% 10*log10

% Frequency range selection
if ~isempty(frange)
    freqRange = f >= frange(1) & f <= frange(2);
else
    freqRange = true(size(f));
end

% Plot
figure('Position',[10 10 1500 700])
hold on
grid on 
[nsignal, ~] = size(xtrue);

for i = 1:nsignal
    [H,f] = tfestimate(xexcit,xtrue, wind, [], nfs, fs);

    psdTrue_dB =  mag2db(abs(H));

    if isdb
        psdTrue_dB =  mag2db(abs(H));

        ylabel("Amplitude [dB]", "Interpreter", "latex")
    else
        psdTrue_dB =  abs(H);

        ylabel("Amplitude [mm]", "Interpreter", "latex")
    end

    plot(f(freqRange), psdTrue_dB(freqRange),  'o', 'MarkerSize', 5, 'MarkerFaceColor','k')
end

% Plot mean of distribution
plot(f(freqRange), yMean_dB(freqRange), 'r', 'LineWidth', 1)

% Patch for confidence interval
patch([f(freqRange); flipud(f(freqRange))], ...
    [yHigh_dB(freqRange); flipud(yLow_dB(freqRange))], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel("Frequency [Hz]", "Interpreter", "latex")
set(gca,'TickLabelInterpreter','latex','fontsize',16, "YScale", scale)
hold off

% Legend
if ~isempty(leg)
    legend(leg, "Interpreter", "latex", "Location", "northwest")
end
set(gca, 'Box', 'on')
if snapShot
    wn = mean(f(freqRange));
    frange            = f>ratio* wn& f<(2-ratio)*wn;
    
    fsnpsht           = f(frange);
    yMean_dB_snpsht   =  (yMean_dB(frange));
    yHigh_dB_snpsht   =  (yHigh_dB(frange));
    yLow_dB_snpsht    =  (yLow_dB (frange));
    psdTrue_dB_snpsht =  (psdTrue_dB(frange) );

    insetPos = [0.65 0.65 0.25 0.26];
    axesInset = axes('Position',insetPos);

    plot(axesInset , fsnpsht, psdTrue_dB_snpsht,   'o', 'MarkerSize', 5, 'MarkerFaceColor','k')
    hold on 
    plot(axesInset, fsnpsht, yMean_dB_snpsht, 'r', 'LineWidth', 1)

    patch(axesInset, [fsnpsht; flipud(fsnpsht)], ...
        [yHigh_dB_snpsht; flipud(yLow_dB_snpsht)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off
    if sum(frange)>1

        xlim(axesInset, [min(fsnpsht) max(fsnpsht)]) 
    end
    grid on
    set(axesInset,'Box','on',"YScale", scale)
    


end

if savePlot
    parent_dir = strcat(pwd, '\Results\',foldname);
    folder_name = fullfile(parent_dir, foldnameSpecial);

    full_path = folder_name;

    if ~isfolder(full_path)

        mkdir(full_path)
    end
    unique_id = string(datetime('now','TimeZone','local','Format','HHmmss'));
    fileName = strcat( "Plot-", "FRF", "-", num2str(unique_id), ".", format);



    folderpath = fullfile(full_path, fileName);
    switch format
        case "pdf"
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperOrientation', 'landscape');
            set(gcf, 'PaperSize',11.7*[1 9/16 ]);
            print(gcf, '-dpdf', '-fillpage', folderpath,  '-vector');
        otherwise
            saveas(gcf, folderpath)
    end
end

end
