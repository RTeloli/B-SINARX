function P = build_embedding_matrix(source_keys, target_keys)
    n_target = length(target_keys);
    n_source = length(source_keys);
    P = zeros(n_target, n_source);
    
    for i = 1:n_target
        key = target_keys{i};
        idx = find(strcmp(key, source_keys));
        if ~isempty(idx)
            P(i, idx) = 1;
        end
    end
end