function f = norm_matrix_rank(bits, fNumber)
%norm_matrix_rank gets the rank of norm matrix for given bits and fermion
%number.
%   return value is in the format [r1, r2, ..., rN], where N = bits

    f = zeros(bits, 1);
    
    for i = 1 : bits
        f(i) = rank(norm_matrix(bits, fNumber, i));
    end
end

