function f = ham4(N)
f=[
    8 24/N 24/N 32/N 0 0 0 4*1i 0 0;
    18/N 6 0 0 32/N 16/N 0 2/N*1i 6*1i 0;
    6/N 0 6 0 0 16/N 0 -2/N*1i 0 2*1i;
    16/N 0 0 8 8/N 8/N 0 -4/N*1i 0 0;
    0 12/N 0 4/N 4 0 24/N 0 0 0;
    0 6/N 18/N 4/N 0 4 24/N 0 -6/N*1i -2/N*1i;
    0 0 0 0 4/N 4/N 0 0 0 0;
    -4*1i -8/N*1i 0 -8/N*1i 0 0 0 0 24/N 8/N;
    -2/N*1i -2*1i 0 0 -8/N*1i 0 0 2/N -6 0;
    -2/N*1i 0 -6*1i 0 0 -8/N*1i 0 2/N 0 -6];
end
