function f=lowest_energies(xi, bits, n, statenumber, mode)
    if nargin < 5
        mode = 'sr';
    end
    
    [Vec,D] = eigs(ham(bits, n) + xi * deltaHam(bits, n), statenumber, mode);
    state_norm = real(diag(ctranspose(Vec)*metric(bits, n)*Vec));
    f = sortrows(horzcat(real(diag(D)),imag(diag(D)), state_norm), 1);
end