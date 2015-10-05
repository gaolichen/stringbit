function f = ground_vanish_N(xi, maxBit, maxN)
    if nargin < 3
        maxN = maxBit;
    end
    
    ret = zeros(maxBit - 2, 5);
    for bit = 3 : maxBit
        ret(bit - 2, 1) = bit;
        for n = min(bit, maxN) : -0.5 : 0.5
            try
                state = lowest_energies(xi, bit, n, 2);
                if state(1, 3) < 1e-8 || abs(state(1, 2)) > 1e-8
                    ret(bit - 2, 2) = n;
                    ret(bit - 2, 3) = state(1, 1);
                    ret(bit - 2, 4) = state(1, 2);
                    ret(bit - 2, 5) = state(1, 3);
                    break;
                end
            catch ME
                fprintf('%s: %s\n', ME.identifier, ME.message);
                fprintf('bit=%d n=%f\n', bit, n);
            end
        end
        
        if ret(bit - 2, 2) == min(bit, maxN)
            for n = min(bit, maxN) + 0.5 : 0.5 : min(bit, maxN) + 5
                try
                    state = lowest_energies(xi, bit, n, 2);
                    if state(1, 3) < 1e-8 || abs(state(1, 2)) > 1e-8
                        ret(bit - 2, 2) = n;
                        ret(bit - 2, 3) = state(1, 1);
                        ret(bit - 2, 4) = state(1, 2);
                        ret(bit - 2, 5) = state(1, 3);
                        break;
                    end
                catch ME
                    fprintf('%s: %s\n', ME.identifier, ME.message);
                    fprintf('bit=%d n=%f\n', bit, n);
                end
            end
        end
    end
    
    f = ret;
end