function mainString             = printFunctions(col, l, Xi)
n = size(Xi, 2);
superscript_map = ['⁰' '¹' '²' '³' '⁴' '⁵' '⁶' '⁷' '⁸' '⁹'];
subscript_map = ['₀' '₁' '₂' '₃' '₄' '₅' '₆' '₇' '₈' '₉'];
for kk = 1:n
    sta = [];
    eta = Xi(: ,kk);
    if all(eta==0)
        sta = "";
    else
        for j= 1:size(l, 1)

            if eta(j)~=0
                if abs(eta(j)) > 100
                    scie = sprintf("%.2e", eta(j));
                    sta = strcat(sta,scie);
                else
                    sta = strcat(sta,num2str(eta(j)));
                end
                for i = 1:length(col)
                    if l(j, i)~=0
                        if l(j, i)== 1
                            ss = sprintf("<strong>%s</strong>" ,    col{i});
                            sta = strcat(sta,ss);
                        else

                            ss = sprintf("<strong>%s%s</strong>" ,    col{i}, superscript_map(l(j, i)+1) );
                            sta = strcat(sta,ss);
                        end
                    end
                end
                if j~=size(l, 1)
                    sta = strcat(sta, " + ");
                end

            end
            if j==size(l, 1)
                a = char(sta);
                if strcmp(a(end-1),'+')
                    a(end-2:end) = [];
                    sta = string(a);
                end

            end

        end
    end
    mainString{kk} = sta;
end
end
