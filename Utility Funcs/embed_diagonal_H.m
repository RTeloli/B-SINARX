function H_M = embed_diagonal_H(H_diag, m1_keys, M_keys)
    % H_diag: vector for the diagonal of H in m1 space
    % m1_keys: keys in m1 (cell array)
    % M_keys: keys in unified space M (cell array)

    nM = length(M_keys);
    H_diag_M = ones(nM, 1);  % default to 1 (multiplicative identity)

    % Create key-to-value mapping
    h_map = containers.Map(m1_keys, H_diag);

    for i = 1:nM
        key = M_keys{i};
        if isKey(h_map, key)
            H_diag_M(i) = h_map(key);
        end
    end

    H_M = diag(H_diag_M);
 end
