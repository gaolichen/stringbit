function [e, f] = get_state(bits, n, index)
    [V,D] = eigs(ham(bits, n), index, 'sr');
    w = V(:, end);
    e = D(index, index);
    f = [];
    for i = 1 : size(w, 1)
        if abs(w(i)) > 1e-8
            f = cat(1, f, [i, w(i)]);
        end
    end
end