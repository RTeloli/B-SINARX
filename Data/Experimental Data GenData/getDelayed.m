function X = getDelayed(dx, delay)

[n ,m] = size(dx);
X = zeros(n-delay,delay*m);
l = 1;
for i=1:delay
    for j = 1:m
        X(:, l) = dx(i:end-delay+(i-1),j );
        l = 1+l;
        
    end
end
end