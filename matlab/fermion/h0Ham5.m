function f = h0Ham5(N)
f=[
    10 16/N 16/N 24/N 24/N 0 0 0 0 0 0 0 4*1i 2*1i 0 0 0 0 0 0 0;
    8/N 8 0 0 0 24/N 12/N 16/N 0 0 0 0 0 4/N*1i 4*1i 0 0 0 0 0 0;
    2/N 0 8 0 0 0 12/N 0 16/N 0 0 0 -2/N*1i 0 0 2*1i 0 0 0 0 0;
    12/N 0 0 10 0 4/N 0 8/N 16/N 0 0 0 -2/N*1i -2/N*1i 0 0 6*1i 0 0 0 0;
    8/N 0 0 0 10 0 4/N 8/N 0 0 0 0 0 -4/N*1i 0 0 0 2*1i 0 0 0;
    0 6/N 0 0 0 6 0 0 0 24/N 8/N 0 0 0 2/N*1i 0 0 0 6*1i 0 0;
    0 2/N 8/N 0 0 0 6 0 0 0 16/N 0 0 0 -2/N*1i 0 0 0 0 2*1i 0;
    0 8/N 0 4/N 6/N 0 0 8 0 12/N 8/N 0 0 0 -4/N*1i 0 0 -2/N*1i 0 0 0;
    0 0 8/N 2/N 0 0 0 0 8 0 4/N 0 0 0 0 -2/N*1i -6/N*1i 0 0 0 0;
    0 0 0 0 0 4/N 0 0 0 4 0 16/N 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 2/N 6/N 0 0 0 4 24/N 0 0 0 0 0 0 -6/N*1i -2/N*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -4*1i -8/N*1i 0 -8/N*1i -12/N*1i 0 0 0 0 0 0 0 -2 4 8/N 4/N 0 0 0 0 10*1i;
    -2*1i -4/N*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 4 2 8/N 4/N 24/N 8/N 0 0 0;
    -4/N*1i -4*1i 0 0 0 -16/N*1i 0 -8/N*1i 0 0 0 0 4/N -4/N 0 0 0 0 24/N 4/N 0;
    -2/N*1i 0 -8*1i 0 0 0 -12/N*1i 0 -16/N*1i 0 0 0 2/N -2/N 0 0 0 0 0 8/N -10/N*1i;
    -2/N*1i 0 0 -2*1i 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 -2 0 4/N 0 -10/N*1i;
    -4/N*1i 0 0 0 -6*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 -2 0 4/N 0;
    0 -2/N*1i 0 0 0 -2*1i 0 0 0 -12/N*1i 0 0 0 0 -2/N 0 0 0 -6 0 0;
    0 -2/N*1i -8/N*1i 0 0 0 -6*1i 0 0 0 -16/N*1i 0 0 0 -2/N -4/N 0 0 0 -6 0;
    0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 -4/N*1i 0 0 -4/N*1i 0 0 -10];
end