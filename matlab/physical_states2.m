function f=physical_states2(xi, bits, n, statenumber, mode, tol)
%physical_states2: Returns statenumber physical states for give xi, bitnumber, n
% If the number of states less than statenumber, return all the states.

    if nargin < 5
        mode = 'sr';
    end
    if nargin < 6
        tol = 10^(-6);
    end
    
    fprintf('xi=%d, bits=%d, n=%d\n', xi, bits, n);
    
    opts.tol = tol;
    tic;
    B = metric(bits, n);
    A = B * (-ham(bits, n) + xi * deltaHam(bits, n));
    %fprintf('a\n');
    idx = licols(B);
    B = B(idx(:), idx(:));
    A = A(idx(:), idx(:));
    %fprintf('b\n');
    if statenumber > size(A, 2)
        statenumber = size(A, 2);
    end
    %disp(A);
    %disp(B);
    
    [Vec,D] = eigs(A, B, statenumber, mode, opts);
    f = diag(D);
end

