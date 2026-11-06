function [contained, stringExpression, containedArray] = getStringInvolvedMain(vars, stringExp)
contained = '';

s = char(stringExp);
s = erase(s, '<strong>');
s = erase(s, '</strong>');
s = regexprep(s, '<[^>]*>', '');            % any leftover HTML tags

if any(double(s) > 128)
    s = nonAscIIReplacement(s);             % z1² -> z1^2
end

% --- implicit multiplication ---
s = regexprep(s, '(\d)([a-df-zA-DF-Z])', '$1*$2');    % v1z1sign -> v1*z1*sign
s = regexprep(s, '(\d)([eE])(?![+\-]?\d)', '$1*$2');  % keep 6.65e-02 intact
s = regexprep(s, '(\))([a-zA-Z(])', '$1*$2');         % sign(z1)v1 -> sign(z1)*v1

% |x| -> abs(x)
s = regexprep(s, '\|\s*([^|]+?)\s*\|', '*abs($1)');

stringExpression = s;

counter = 0;
containedArray = zeros(size(vars));
for i = 1:length(vars)
    pat = ['(?<![A-Za-z_])' regexptranslate('escape', char(vars{i})) '(?!\w)'];
    if ~isempty(regexp(s, pat, 'once'))
        containedArray(i) = 1;
        if counter == 0
            contained = vars{i};
        else
            contained = strcat(contained, ",", vars{i});
        end
        counter = counter + 1;
    end
end

    function z = nonAscIIReplacement(z)
        superscript_map = ['⁰' '¹' '²' '³' '⁴' '⁵' '⁶' '⁷' '⁸' '⁹'];
        dm = '0123456789';
        r = unique(z(double(z) > 128));
        for k = 1:numel(r)
            idx = find(superscript_map == r(k), 1);
            if ~isempty(idx)
                z = strrep(z, r(k), ['^' dm(idx)]);
            end
        end
    end
end