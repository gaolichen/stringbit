function f = ground_vanish_N(maxBit)
    ret = zeros(maxBit - 2, 2);
    for bit = 3 : maxBit
        ret(bit - 2, 1) = bit;
        for n = bit : -0.5 : 0.5
            state = lowest_energies(bit, n, 1);
            if state(1, 3) < 1e-8 || abs(state(1, 2)) > 1e-8
                ret(bit - 2, 2) = n;
                break;
            end
        end
    end
    
    f = ret;
end