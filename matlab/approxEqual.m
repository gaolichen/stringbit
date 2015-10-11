function f = approxEqual(a, b)
    f = abs(a - b) < 1e-3;
end