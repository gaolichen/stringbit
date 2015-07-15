function f=lowest_energies2(bits, n, statenumber)
    [Vec,D] = eigs(pham(bits, n), statenumber, 'sr');
    state_norm = real(diag(ctranspose(Vec)*metric(bits, n)*Vec));
    f = sortrows(horzcat(real(diag(D)),imag(diag(D)), state_norm), 1);
end