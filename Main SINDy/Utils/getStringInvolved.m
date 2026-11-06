function [contained, stringExpression] = getStringInvolved(vars, stringExp)

contained = '';
stirngexp = char(stringExp);
stringExpression = stirngexp;
counter = 0; 
for i = 1:length(vars)
    array  = contains(stirngexp,vars{i});
    if array
        if counter ==0
            contained = vars{i};
            counter = counter +1; 
        else 
            contained = strcat(contained, ",",vars{i});
            counter = counter +1; 
        end

    stringExpression= strrep(stringExpression, string( vars{i}), strcat("*", string( vars{i})));
    end
end
stringExpression = erase(stringExpression, "<strong>"); 
stringExpression = erase(stringExpression, "</strong>"); 
% end