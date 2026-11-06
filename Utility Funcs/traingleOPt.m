
%
function maxIndex= traingleOPt(f1, f2)

if  checkIfTri(f1, f2) 
 
    p11 = f1(1);
    p21 = f2(1);

    p12 = f1(end);
    p22 = f2(end);
    p1 = [ p11 p21  ];

    p2 = [ p12 p22 ];
    area = zeros(length(f1), 1);
    for i = 1:length(f1)
        p3 = [f1(i) f2(i)];
        area(i) = AreaTriangle(p1, p2, p3);
    end

    % [~, maxIndex] = max(area);
    maxIndex = find(area==max(area));
    maxIndex = maxIndex(end);
    p3 = [f1(maxIndex) f2(maxIndex)];
else
    [~, maxIndex] = max(f1);
end 

end
function area = AreaTriangle(p1, p2, p3)


area = [[p1 1]' [p2 1]' [p3 1]'];
area = 0.5*abs(det(area));
end

function indc = checkIfTri(f1, f2)



indicator1 = (max(f1)-min(f1))/max(f1);
indicator2 = (max(f2)-min(f2))/max(f2);

thres =1;
indc = abs(indicator1*100)>thres  && abs(indicator2*100)>thres ;
end