function [M , P, hUpdated]= createSimilarityMap(B, S, xis, xib, H)
m1 = containers.Map(B, xib);
m2 = containers.Map( S , xis);



all_keys = unique([B; S]);

M_keys = all_keys;


M_values = zeros(1, length(all_keys));
M_source = cell(1, length(all_keys));


for i = 1:length(all_keys)
    key = all_keys{i};
    in_m1 = isKey(m1, key);
    in_m2 = isKey(m2, key);

    val = 0;

    if in_m1
        val = val + m1(key);
    end
    if in_m2
        val = val + m2(key);
    end

    M_values(i) = val;

    if in_m1 && in_m2
        M_source{i} = 'both';
    elseif in_m1
        M_source{i} = 'm1';
    else
        M_source{i} = 'm2';
    end
end
 disp(table(M_keys, M_values', M_source', ...
    'VariableNames', {'Key', 'Value', 'Source'}));
M = containers.Map(M_keys, M_values);



nM = length(M_keys);
 

n1     = length(B);
P      = zeros(nM, n1);
m1_map = containers.Map(B, 1:n1);

for i = 1:nM
    key = M_keys{i};
    if isKey(m1_map, key)
        j = m1_map(key);
        P(i, j) = 1;
    end
end 

H_diag = diag(H); 
hUpdated = embed_diagonal_H(H_diag, B, M_keys);
end

 