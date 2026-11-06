function [sol ,l,   string, mainString]                     = lib_1(x, p, vars)

[n, ~] = size(x) ;

l = find_combinations_efficient(n, p, p) ;
sol = [];
string =  vars;
for i = 1:size(l, 1)
    mul = x'.^l(i, :);
    if n>1
        mul = prod(mul')';
    end
    if isempty(sol)
        mul = ones(size(mul, 1),1);
    end
    sol = [sol mul];
end


L = l;
superscript_map = ['⁰' '¹' '²' '³' '⁴' '⁵' '⁶' '⁷' '⁸' '⁹'];

sta = {};
for i = 1:size( L ,1)
    s = "";
    for j = 1:size( L ,2)
        if L(i, j)~=0
            if L(i, j)==1

                s = strcat(s, string{j});

            else
                s = strcat(s, string{j},superscript_map(L(i,j)+1));
            end
        end

    end
    sta{i} = s;
end
mainString = [sta{:}]';
end
