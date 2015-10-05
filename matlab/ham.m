function f = ham(bits, n)
    if bits == 3
        f = h0Ham3(n);
    elseif bits == 4
        f = h0Ham4(n);
    elseif bits == 5
        f = h0Ham5(n);
    elseif bits == 6
        f = h0Ham6(n);
    elseif bits == 7
        f = h0Ham7(n);
    elseif bits == 8
        f = h0Ham8(n);
    elseif bits == 9
        f = h0Ham9(n);
    elseif bits == 10
        f = h0Ham10(n);
    elseif bits == 11
        f = h0Ham11(n);
    else
        f = 0;
    end
end