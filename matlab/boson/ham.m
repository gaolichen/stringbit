function f = ham(bits, n)
    if bits == 3
        f = ham3(n);
    elseif bits == 4
        f = ham4(n);
    elseif bits == 5
        f = ham5(n);
    elseif bits == 6
        f = ham6(n);
    elseif bits == 7
        f = ham7(n);
    elseif bits == 8
        f = ham8(n);
    elseif bits == 9
        f = ham9(n);
    elseif bits == 10
        f = ham10(n);
    elseif bits == 11
        f = ham11(n);
    else
        f = 0;
    end
end