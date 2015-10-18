function f=physical_states(xi, bits, n, statenumber, mode, tol)
    if nargin < 5
        mode = 'sr';
    end
    if nargin < 6
        tol = 10^(-6);
    end
    
    opts.tol = tol;
    tic;
    B = metric(bits, n);
    disp(toc);
    A = B * (ham(bits, n) + xi * deltaHam(bits, n));
    disp(toc);
    idx = licols(B);
    disp(toc);
    B = B(idx(:), idx(:));
    A = A(idx(:), idx(:));
    disp(toc);
    %disp(A);
    %disp(B);
    disp(size(A));
    disp(size(B));
    [Vec,D] = eigs(A, B, statenumber, mode, opts);
    f = diag(D);
    disp(toc);
end
