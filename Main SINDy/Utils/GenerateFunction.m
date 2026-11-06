function [mainString, evalString] = GenerateFunction(funcs, Xi)
[m, n] = size(Xi);
mainString = cell(1, n);
evalString = cell(1, n);

for kk = 1:n
    eta = Xi(:, kk);
    dispTerms = {};
    evalTerms = {};
    signs     = [];

    if all(eta == 0)
        mainString{kk} = "";
        evalString{kk} = "0";
        continue;
    end

    for j = 1:m
        if eta(j) == 0
            continue;
        end

        % coefficient magnitude (sign handled separately when joining)
        c = abs(eta(j));
        if c > 100 || c < 1e-3
            coefStr = sprintf('%.4e', c);
        else
            coefStr = num2str(c);
        end

        % --- sanitize term label ---
        term = char(funcs{j});
        term = regexprep(term, '^\*+', '');          % stray leading '*'
        term = regexprep(term, '\*+$', '');          % stray trailing '*'
        term = regexprep(term, '(?<=\d)(?=[a-zA-Z])', '*');  % v1z1 -> v1*z1, z1sign -> z1*sign
        term = regexprep(term, '(?<=\))(?=[a-zA-Z(])', '*'); % )( or )x  -> )*( , )*x
        term = regexprep(term, '\*{2,}', '*');       % collapse accidental '**'

        signs(end+1)     = sign(eta(j));                                   %#ok<AGROW>
        dispTerms{end+1} = sprintf('%s*<strong>%s</strong>', coefStr, term); %#ok<AGROW>
        evalTerms{end+1} = sprintf('%s*%s', coefStr, term);                 %#ok<AGROW>
    end

    % --- join terms with correct signs ---
    if signs(1) < 0
        dispStr = strcat("-", dispTerms{1});
        evalStr = strcat("-", evalTerms{1});
    else
        dispStr = string(dispTerms{1});
        evalStr = string(evalTerms{1});
    end

    for t = 2:numel(dispTerms)
        if signs(t) < 0
            op = " - ";
        else
            op = " + ";
        end
        dispStr = strcat(dispStr, op, dispTerms{t});
        evalStr = strcat(evalStr, op, evalTerms{t});
    end

    mainString{kk} = dispStr;
    evalString{kk} = evalStr;
end
end