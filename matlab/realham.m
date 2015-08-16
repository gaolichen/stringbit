function f = realham(bits, n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    H = ham(bits, n);
    G = metric(bits, n);
    %S = orthnormal(G, true);
    S = orthnormal(G);
    f = S * G * H * ctranspose(S);
end

