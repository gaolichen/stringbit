function f = deltaHam(bits, n)
    if bits == 3
        f = deltaHam3(n);
    elseif bits == 4
        f = deltaHam4(n);
    elseif bits == 5
        f = deltaHam5(n);
    elseif bits == 6
        f = deltaHam6(n);
    elseif bits == 7
        f = deltaHam7(n);
    elseif bits == 8
        f = deltaHam8(n);
    elseif bits == 9
        f = deltaHam9(n);
    elseif bits == 10
        f = deltaHam10(n);
    elseif bits == 11
        f = deltaHam11(n);
    else
        f = 0;
    end
end