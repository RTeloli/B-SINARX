 
function [sol,  vars, mainStirng,l] = libSignPoly(x, p, vars)
        n = size(x, 1);
        sign_functions = sign(x);
        fullform = [x; sign_functions];
        string = getSignString(vars);
        [sol, l] = lib2(fullform, p, [vars vars]);
        L = l;
        superscript_map = ['⁰' '¹' '²' '³' '⁴' '⁵' '⁶' '⁷' '⁸' '⁹'];
        sta = {};
        for i = 1:size(L, 1)
            s = "";
            for j = 1:size(L, 2)
                if L(i, j) ~= 0
                    if L(i, j) == 1
                        s = strcat(s, string{j});
                    else
                        s = strcat(s, string{j}, superscript_map(L(i,j)+1));
                    end
                end
            end
            sta{i} = s;
        end
        clc
        mainStirng = [sta{:}]';
    function col = getSignString(col)
        sign_col = {};
        for i = 1:length(col)
            sign_col{i} = strcat('sign(', col{i}, ')');
        end
        col = [col sign_col];
    end
end