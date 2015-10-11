function f=lowest_energies2(xi, bits, n, statenumber, mode, tol)
    if nargin < 5
        mode = 'sr';
    end
    if nargin < 6
        tol = 10^(-6);
    end
    
    opts.tol = tol;
    [Vec,D] = eigs(-ham(bits, n) + xi * deltaHam(bits, n), statenumber, mode, opts);
    state_norm = real(diag(ctranspose(Vec)*metric(bits, n)*Vec));
    f = sortrows(horzcat(real(diag(D)),imag(diag(D)), state_norm), 1);
end