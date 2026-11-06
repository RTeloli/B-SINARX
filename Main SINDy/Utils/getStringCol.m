function col                    = getStringCol(n,vars)

if strcmp(vars,"")
    vars = cellstr(repmat("x", n, 1));
    for i = 1:n
        col{i} = strcat(vars{i},'_', num2str(i));
    end
else
    for i = 1:n
        col{i} = strcat(vars{i},'', "");
    end

end 
end