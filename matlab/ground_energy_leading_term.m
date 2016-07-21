function f = ground_energy_leading_term(xi, bits, n, tol)
%ground_energy_leading_term If the ground energy is const + a/N^2 + higher
%order term, then the function returns a.
%   Detailed explanation goes here
    if nargin < 4
        tol = 10^(-6);
    end

a = lowest_energies(xi, bits, n, 3, 'sr', tol);
b = lowest_energies(xi, bits, 2*n, 3, 'sr', tol);
c = (a - b) * 4 * n^2 / 3;
f = c(1,1);

end