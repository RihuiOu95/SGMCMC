function [result] = identify_brust_region(y, cutoff, gap, left)
    if nargin <= 3
        left = false;
    end
    if left == false
        indices = find(y>cutoff);
    else
        indices = find(y<cutoff);
    end
    result = [];
    for i = 2 : length(indices)
        if (indices(i)>indices(i-1)+gap)
            result = [result, indices(i)];
        end
    end
end