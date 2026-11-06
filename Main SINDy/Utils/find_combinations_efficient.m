function valid_combinations           = find_combinations_efficient(n, m, p)
% n: number of spots
% m: range of numbers (0 to m)
% p: maximum sum

% Range of numbers
range = 0:m;

% Initialize results
valid_combinations = [];

% Recursive helper function
    function generate_combination(current_comb, remaining_sum, remaining_spots)
        if remaining_spots == 0
            % Accept the combination if the remaining_sum is non-negative
            if remaining_sum >= 0
                valid_combinations = [valid_combinations; current_comb];
            end
            return;
        end
        for i = range
            if remaining_sum - i >= 0
                generate_combination([current_comb, i], remaining_sum - i, remaining_spots - 1);
            end
        end
    end

% Start recursion
generate_combination([], p, n);
end
