function f = orthnomal(A)
    B = orth(A);
    C = ctranspose(B);
    D = C * A * ctranspose(C);
    for i = 1 : size(A, 1)
        C(i,:) = C(i,:) / sqrt(abs(D(i, i)));
    end
    
    f = C;
end