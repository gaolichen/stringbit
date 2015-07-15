function f = metric(bits, n)
    if bits == 3
        f = norm3(n);
    elseif bits == 4
        f = norm4(n);
    elseif bits == 5
        f = norm5(n);
    elseif bits == 6
        f = norm6(n);
    elseif bits == 7
        f = norm7(n);
    elseif bits == 8
        f = norm8(n);
    elseif bits == 9
        f = norm9(n);
    elseif bits == 10
        f = norm10(n);
    elseif bits == 11
        f = norm11(n);
    else
        f = 0;
    end
end