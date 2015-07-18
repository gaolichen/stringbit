function [e, f0, f1] = get_ground_state(bits, n)
    [V,D] = eigs(ham(bits, Inf), 1, 'sr');
    w0 = V(:, end);
    d0 = D(1, 1);
    
    [V1,D1] = eigs(ham(bits, n), 1, 'sr');
    w1 = V1(:, end);
    d1 = D1(1, 1);
    
    [V2,D2] = eigs(ham(bits, 2 * n), 1, 'sr');
    w2 = V2(:, end);
    d2 = D2(1, 1);
    
    e0 = d0;
    e1 = n * (4 * d2 - d1 - 3 * d0);
    
    e = [e0 e1];
    f0 = [];
    f1 = [];
    
    for i = 1 : size(w0, 1)
        if abs(w0(i)) > 1e-8
            f0 = cat(1, f0, [i, w0(i), abs(w0(i))]);
        else
            res = n * (4 * w2(i) - w1(i) - 3 * w0(i));
            if abs(res) > 0.01
                f1 = cat(1, f1, [i, res, abs(res)]);
            end
        end
    end
    
    f0 = flipud(sortrows(f0, 3));
    f1 = flipud(sortrows(f1, 3));
end