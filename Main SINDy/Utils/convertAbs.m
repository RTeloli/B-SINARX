function [convertedAbs]  = convertAbs(expres)
convertedAbs = regexprep(expres, 'abs\((.*?)\)', '|$1|');
end