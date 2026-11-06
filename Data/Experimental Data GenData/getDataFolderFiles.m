
function [Z, U, V, X , A, t, params, dz ] = getDataFolderFiles(  varargin)

[U, V, X , A, t ] = getDataFull( varargin{:} );
 
zeta = 0.00446;     %
mu_xi   =  0.444/100 ;
co_xi   = 0.57/100;
mu_fn   =  18.80;
wn      = 2*pi*mu_fn;
alpha = 5.000292979243245e+02;
k = wn^2 - alpha;
m = 1-alpha/wn^2;
U = 49.5*U; 

Z  = U  -  (A + zeta*2*wn*V +wn^2*X)  ;

U = U(:); 
V = V(:); 
X = X(:); 
A = A(:); 
Z = Z(:); 
dz = gradient(Z, t);
params.m = 1; 
params.c = 2*mu_xi*wn;   
params.k   = wn^2;  
params.kk   = k;  
params.alpha  = alpha;  

end
function  [U_filter , v_filter, Disp , Acc, t  ] = getDataFull( varargin)

SignalType      = "Chirp";
index           = 2;
excitationLevel = "30"; 
nExperNeeded = 1;
mainFolder = "Data\Experimental Data GenData\Data Folder\";
selection = 0; 
for i  = 1:length(varargin)
    if isstring(varargin{i})
        switch varargin{i}
            case "n_experiments"
                nExperNeeded    = varargin{i+1};
            case "amplitude"
                excitationLevel = varargin{i+1};
            case "cuttoffFreq"
                fc              = varargin{i+1};
            case "index"
                index           = varargin{i+1};
            case "signal"  
                SignalType      = varargin{i+1};
            case "selection"
                selection       = varargin{i+1};

        end
    end
end

try
    switch SignalType.lower()
        case "chirp"
            reconName         = strcat(SignalType, "0p",excitationLevel);
            loadFOlder        = strcat(mainFolder, SignalType);
            loadFOlder        = fullfile(pwd, loadFOlder);

            listName0         = ls(loadFOlder);
            listName0(1:2, :) = [];
            selectedSub       = listName0(index,:);
            subfolder         = strcat(loadFOlder, "\", selectedSub);
            listName          = ls(subfolder);
            listName(1:2, :)  = [];
            fileFolder        = strcat(subfolder, "\",reconName);

            checkList         = string(ls(fileFolder));
            searchedList      = strfind(checkList,"Time" ) ;

            isOne             = cellfun(@(x) ~isempty(x), searchedList);
            isOneList         = checkList(isOne, :);
            data              = struct();

            for i = 1:size(checkList(isOne, :), 1)
                fileName = isOneList(i);
                fullName= fullfile(fileFolder,fileName );
                temp = load(fullName);
                fieldN = fieldnames(temp);
                for j = 1:length(fieldN)
                    field = fieldN{j};
                    if contains(field, "Time")
                        field_new= replace(field,"Time", "n");
                    else
                        field_new = field;
                    end
                    data.(field_new) = temp.(field);
                end

            end

            v  = data.n_Laser__X.y_values.values;
            u  = data.n_Shaker__X.y_values.values;
            t0 = data.n_Laser__X.x_values.start_value;
            h  = data.n_Laser__X.x_values.increment;
            N  = data.n_Laser__X.x_values.number_of_values;

            Ntotal =  length(u);
            Neach = 16384;
            r = floor(Ntotal/Neach);
            U = reshape(u(1:r*Neach), [Neach r]);
            V = reshape(v(1:r*Neach), [Neach r]);
            t = t0:h:(Neach-1)*h+t0;




        case "rand"
            reconName = strcat(SignalType, "0p",excitationLevel);
            loadFOlder= strcat(mainFolder, SignalType);
            loadFOlder = fullfile(pwd, loadFOlder);

            subfolder = strcat(loadFOlder, "\", reconName);
            checkList = string(ls(subfolder));
            searchedList = strfind(checkList,"Time" ) ;
            isOne = cellfun(@(x) ~isempty(x), searchedList);
            isOneList = checkList(isOne, :);

            data = struct();

            for i = 1:size(checkList(isOne, :), 1)
                fileName = isOneList(i);
                fullName= fullfile(subfolder,fileName );
                temp = load(fullName);
                fieldN = fieldnames(temp);
                for j = 1:length(fieldN)
                    field = fieldN{j};
                    if contains(field, "Time")
                        field= replace(field,"Time", "n");
                    end
                    data.(field) = temp.(field);
                end

            end

            v = data.n_Laser__X.y_values.values;
            u  = data.n_Shaker__X.y_values.values;
            t0 = data.n_Laser__X.x_values.start_value;
            h  = data.n_Laser__X.x_values.increment;
            N  = data.n_Laser__X.x_values.number_of_values;

            Ntotal =  length(u);
            rr = (abs(u-0.2)<1e-4);
            Neach = 16384;
            r = floor(Ntotal/Neach);
            U = reshape(u(1:r*Neach), [Neach r]);
            V = reshape(v(1:r*Neach), [Neach r]);
            t = t0:h:(Neach-1)*h+t0;
        case "sine"
            reconName = strcat(SignalType, "_0p",excitationLevel, "Hz");
            loadFOlder= strcat(mainFolder,  SignalType);
            loadFOlder = fullfile(pwd, loadFOlder);
            listName0 = ls(loadFOlder);
            listName0(1:2, :) = [];
            reconName = listName0(index,:);

            subfolder = strcat(loadFOlder, "\", reconName);
            checkList = string(ls(subfolder));
            searchedList = strfind(checkList,"Time" ) ;
            isOne = cellfun(@(x) ~isempty(x), searchedList);
            isOneList = checkList(isOne, :);

            data = struct();

            for i = 1:size(checkList(isOne, :), 1)
                fileName = isOneList(i);
                fullName= fullfile(subfolder,fileName );
                temp = load(fullName);
                fieldN = fieldnames(temp);
                for j = 1:length(fieldN)
                    field = fieldN{j};
                    if contains(field, "Time")
                        field= replace(field,"Time", "n");
                    end
                    data.(field) = temp.(field);
                end

            end
            v  = data.n_Laser__X.y_values.values;
            u  = data.n_Shaker__X.y_values.values;
            t0 = data.n_Laser__X.x_values.start_value;
            h  = data.n_Laser__X.x_values.increment;
            N  = data.n_Laser__X.x_values.number_of_values;

            Ntotal =  length(u);
            g2 =  gradient(gradient(gradient(u)));

            [~, peackINdex] = findpeaks(g2, 'MinPeakHeight', mean(g2) + 5*std(g2), 'MinPeakWidth',2);
            c = length(find(peackINdex));

            b =  peackINdex(1:2*floor(c/2));

            d     = unique(diff(b));
            Neach = d;
            t = t0:h:(Neach-1)*h+t0;

            r = floor(Ntotal/Neach);

            U = reshape(u(1:r*Neach), [Neach r]);
            V = reshape(v(1:r*Neach), [Neach r]);

    end


    nEx = nExperNeeded;
    if nExperNeeded>r
        nEx = r;
        fprintf("Maximum experiments allowed is %i\n", r);
    end
    if selection==0

        chosenIndex = randperm(r, nEx); 
    elseif nEx==1
        chosenIndex =selection;
    else
        chosenIndex = selection;
    end
    U = U(:,chosenIndex );
    V = V(:,chosenIndex );
    U_filter = U; %butterWorthFilterHP(U, h, fc);
    v_filter = butterWorthFilterHP(V, h, fc);

    U_filter = U; %butterWorthFilterLP(U_filter, h, 1);
    v_filter = butterWorthFilterLP(v_filter, h, 1);

    x0  = 0;
    Disp = integrateNumVelo(v_filter,t, x0');
    Acc  = zeros(size(v_filter));
    for i = 1:size(v_filter, 2)
        Acc (:, i)   = gradient(v_filter(:, i),t );
    end
    thalf = t(round(length(t)/2));
    t1sth = t(round(length(t)/10));
    tind = t<thalf&t>t1sth;


    U_filter    = U_filter(tind, :);
    v_filter    = v_filter(tind, :);
    Disp        =     Disp(tind, :);
    Acc         =     Acc(tind,  :);
    t           =        t(tind);
catch e
    disp(e.message)
    disp("line "+ e.stack(1).line)
    U_filter = [];
    v_filter = [];
    Disp     = [];
    Acc      = [];
end
end