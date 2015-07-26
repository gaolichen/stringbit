function f = ham4(N)
f=[
    8 36/N 36/N 48/N 0 0 0 4*1i 0 0;
    30/N 6 0 0 48/N 24/N 0 2/N*1i 6*1i 0;
    10/N 0 6 0 0 24/N 0 -2/N*1i 0 2*1i;
    24/N 0 0 8 12/N 12/N 0 -4/N*1i 0 0;
    0 20/N 0 8/N 4 0 36/N 0 0 0;
    0 10/N 30/N 8/N 0 4 36/N 0 -6/N*1i -2/N*1i;
    0 0 0 0 8/N 8/N 0 0 0 0;
    -4*1i -8/N*1i 0 -8/N*1i 0 0 0 0 36/N 12/N;
    -2/N*1i -2*1i 0 0 -8/N*1i 0 0 6/N -6 0;
    -2/N*1i 0 -6*1i 0 0 -8/N*1i 0 6/N 0 -6];
end
