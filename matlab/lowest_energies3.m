function f = lowest_energies3(bits, n, statenumber)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [Vec,D] = eigs(ham(bits, n), statenumber, 'sr');
    G = metric(bits, n);
    state_norm = real(diag(ctranspose(Vec)*G*Vec));
    S = orthnormal(G);
    %disp(S);
    state_norm2 = real(diag(ctranspose(Vec)*G*ctranspose(S)*S*G*Vec)) - state_norm;
    f = sortrows(horzcat(real(diag(D)),imag(diag(D)), state_norm, state_norm2), 1);
end