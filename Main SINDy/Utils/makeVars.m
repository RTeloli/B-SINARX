function [vars, subFunc, storeIndex] = makeVars(eqivars,expre )
n = length(eqivars); 
genericVars = [char(120:122) char((120:122)-3)]; 
subFunc = expre; 
storIndex = []; 
charList = cellstr(num2cell(genericVars'));
for i = 1:n
    logArr = cellfun(@(s) contains(eqivars{i}, s), charList); 
    vars{i} = charList{logArr};  
    storIndex = [storIndex; strfind(subFunc, eqivars{i})]; 
    subFunc = strrep(subFunc, eqivars{i}, vars{i});

end
end

