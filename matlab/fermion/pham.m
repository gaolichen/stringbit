function f = pham(bits, n)
    if bits == 3
        f = pham3(n);
    elseif bits == 4
        f = pham4(n);
    elseif bits == 5
        f = pham5(n);
    elseif bits == 6
        f = pham6(n);
    elseif bits == 7
        f = pham7(n);
    elseif bits == 8
        f = pham8(n);
    elseif bits == 9
        f = pham9(n);
    elseif bits == 10
        f = pham10(n);
    elseif bits == 11
        f = pham11(n);
    else
        f = 0;
    end
end