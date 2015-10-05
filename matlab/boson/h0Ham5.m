function f = h0Ham5(N)
f=[
    10 16/N 24/N 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0;
    10/N 8 0 24/N 16/N 0 0 0 2/N*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0;
    20/N 0 10 4/N 16/N 0 0 -2/N*1i -2/N*1i 0 0 2*1i 0 0 0 0 0 0 0 0 0;
    0 8/N 0 6 0 24/N 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0;
    0 8/N 6/N 0 8 12/N 0 0 0 -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0;
    0 0 0 6/N 0 4 40/N 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -10*1i -16/N*1i -24/N*1i 0 0 0 0 2 4 12/N 4/N 8/N 0 0 0 0 0 4*1i 0 0 0;
    0 0 0 0 0 0 0 4 -2 4/N 0 16/N 4/N 0 0 0 0 -2*1i 0 0 0;
    -10/N*1i -8*1i 0 -24/N*1i -16/N*1i 0 0 2/N -6/N 0 0 0 0 16/N 4/N 0 0 -2/N*1i 0 0 0;
    0 0 0 0 0 0 0 4/N 0 0 8 0 0 0 12/N 16/N 0 0 4*1i 0 0;
    -10/N*1i 0 -6*1i 0 -16/N*1i 0 0 0 0 0 0 -2 0 4/N 0 4/N 0 -6/N*1i 0 0 0;
    0 0 0 0 0 0 0 0 4/N 0 0 0 10 0 4/N -8/N 0 0 0 6*1i 0;
    0 -8/N*1i 0 -6*1i 0 -24/N*1i 0 0 0 -4/N 0 0 0 -6 0 0 4/N 0 0 0 0;
    0 0 0 0 0 0 0 0 0 4/N 6/N 0 0 0 6 0 16/N 0 2/N*1i 0 6*1i;
    0 0 0 0 0 0 0 0 0 0 8/N 4/N -2/N 0 0 8 4/N 0 -4/N*1i 6/N*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 4/N 0 4 0 0 0 0;
    0 0 0 0 0 0 0 -4*1i 2*1i -8/N*1i 0 -8/N*1i 0 0 0 0 0 -10 4/N -12/N 0;
    0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 -8/N*1i -8/N*1i 0 4/N 0 0 12/N;
    0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 -4/N 0 -2 4/N;
    0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -2*1i 0 -8/N*1i 0 -2/N 0 -6];
end
