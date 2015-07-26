function f = ham7(N)
f=[
    14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 0 0 0 0 2*1i 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -12/N -4/N 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 -10/N -6/N 0 12 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 2/N*1i 0 0 0 -2/N*1i 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 -8/N 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 -10/N -4/N 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 -8/N -12/N -8/N 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 -6/N 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -8/N -4/N 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 -6/N -12/N 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -6/N -8/N 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14*1i -24/N*1i -40/N*1i -48/N*1i 0 0 0 0 0 0 0 0 0 0 0 6 4 0 0 0 4/N -8/N 0 0 -12/N 0 -16/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 6 4 0 0 0 16/N -8/N 4/N 12/N 0 32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 2*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 2 0 0 0 -8/N 24/N 0 24/N 4/N -16/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i -2*1i 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i -12*1i 0 0 -40/N*1i -32/N*1i -36/N*1i 0 0 0 0 0 0 0 0 -14/N -6/N 4/N 4 4 0 0 0 0 0 0 0 0 0 4/N -8/N 0 -12/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 2/N*1i 4*1i 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -12/N -12/N 4 4 0 0 0 0 0 0 0 0 0 0 16/N 4/N 24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 -2/N*1i 0 0 0 0 4*1i 2*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i 0 -10*1i 0 0 -16/N*1i 0 -48/N*1i 0 0 0 0 0 0 0 -4/N 0 -8/N 0 0 0 6 4 0 0 0 0 0 0 0 0 0 0 4/N 0 0 -16/N 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i -2/N*1i -4/N*1i -2/N*1i 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -8/N 4/N 0 0 0 4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 48/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i -2/N*1i -4/N*1i -2/N*1i 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 4*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i 0 0 -8*1i 0 0 -24/N*1i -16/N*1i 0 0 0 0 0 0 0 -4/N -4/N 4/N 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -2/N*1i -2/N*1i -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i 0 0 -6*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 0 -8/N 4/N 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i -2/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -12/N*1i 0 0 -10*1i 0 0 0 -48/N*1i -24/N*1i 0 0 0 0 0 0 0 0 -12/N -6/N 0 -4/N 0 0 0 0 0 2 4 0 0 0 0 0 0 0 0 0 0 4/N -8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -14/N 0 0 -4/N 0 0 0 0 4 -2 0 0 0 0 0 0 0 0 0 0 0 24/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 -10/N 0 0 -4/N 0 0 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 2/N*1i 2/N*1i 6/N*1i 0 0 0 0 0 0 0 4*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -12/N*1i -10/N*1i 0 0 -8*1i 0 0 0 -24/N*1i -48/N*1i 0 0 0 0 0 0 0 -4/N -4/N 0 -10/N -10/N 0 -6/N 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 -8/N*1i 0 -4/N*1i 0 0 0 0 -2/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 -8/N 0 -4/N 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -12/N*1i 0 -8/N*1i 0 0 -12*1i 0 0 -16/N*1i 0 0 0 0 0 0 0 0 -8/N 4/N 0 0 0 0 -12/N 0 -8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 -2/N 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N 0 0 4/N 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 2/N 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i -2/N*1i -2/N*1i -6/N*1i 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 6*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 -10/N*1i 0 0 0 0 -6*1i 0 0 -24/N*1i 0 0 0 0 0 0 0 0 0 0 -8/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 -10/N*1i 0 0 0 -8*1i 0 0 -48/N*1i -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -10/N -10/N 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 -8/N 0 -4/N 0 -4/N 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 -10/N*1i -8/N*1i 0 0 0 -6*1i 0 0 -48/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N 4/N 0 -12/N 0 -6/N 0 0 0 -8/N 0 0 0 -2 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 -6/N 0 0 0 0 0 -4/N 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 -6/N 0 -6/N 0 -4/N 0 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 6/N*1i 0 2/N*1i 0 0 -2/N*1i 0 0 0 6*1i 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N 4/N 0 -8/N -4/N 0 0 0 0 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 2/N 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i -2/N*1i 0 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -8/N*1i 0 0 -6*1i 0 -40/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -12/N 0 -4/N 0 0 0 0 -6 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -6/N 0 -4/N -4/N 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 6*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 2/N -4/N -6/N -8/N 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 6/N*1i 0 -2/N*1i 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N -4/N 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -2*1i -8/N*1i -4/N*1i 0 -8/N*1i 0 0 -12/N*1i 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 2 0 2 0 0 0 0 4/N 0 0 0 -8/N 0 0 0 -4/N 12/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 2*1i 0 -4/N*1i 0 0 -8/N*1i 4/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 2 -2 -2 0 0 0 0 0 0 0 16/N 8/N 4/N -6/N 4/N 18/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 6/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 -8/N*1i 0 0 -16/N*1i 0 0 -12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 2 4 0 0 0 0 4/N 4/N 0 -8/N -4/N 0 -4/N -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 2*1i 0 -4/N*1i 0 0 -8/N*1i -4/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 2 -2 -2 0 0 0 0 0 0 0 16/N 0 -4/N -2/N 12/N -18/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 -6/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 -4/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 4 -2 -6 0 0 0 0 0 0 12/N 16/N 0 -4/N 0 8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 2/N*1i -4*1i 0 0 0 0 0 0 0 0 -16/N*1i -8/N*1i 0 -8/N*1i 0 -12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -10/N -8/N -6/N -8/N 0 -4 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 -12/N 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N 0 8/N 0 0 -4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N -24/N -8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 -8/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N -8/N -2/N -6/N 0 0 -4 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 -8/N*1i 0 0 0 -8/N*1i -16/N*1i -12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -2/N 0 2/N 0 0 0 0 0 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N -8/N -36/N -4/N 0 0 0 0 0 0 0 4*1i 0 0 0 -12/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 -8/N*1i 0 -6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N 2/N 2/N 2/N 0 0 0 0 2 0 4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N 18/N -6/N 0 0 0 0 0 0 -2/N*1i 2*1i 0 0 0 6/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 -8/N*1i 0 -6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N -2/N 6/N -2/N 0 0 0 0 2 4 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18/N 2/N 0 0 0 0 0 0 2/N*1i 2*1i 0 0 0 6/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N 0 2/N 4/N 0 0 0 0 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N 0 0 8/N 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 -4*1i 2*1i 0 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -8/N 4/N -4/N 4/N 0 0 0 0 0 0 0 0 -6 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 -12/N 0 0 0 0 0 -10/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 6/N -4/N 2/N 0 0 0 0 0 0 0 0 0 2 4 0 0 0 0 0 0 0 0 0 0 0 0 0 -24/N 0 0 0 0 0 0 0 10*1i 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 2/N 0 -2/N -4/N 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 24/N 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N -8/N 0 0 0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -10/N*1i 2/N*1i 2/N*1i 0 0 0 0 0 0 -6*1i 0 -8*1i 0 0 0 -8/N*1i 0 -12/N*1i 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -12/N -8/N -12/N 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 -6 0 0 0 0 0 0 0 0 4/N 0 4/N 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 4/N 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 0 0 0 -4*1i 2*1i 0 0 0 0 0 0 0 0 0 -24/N*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -12/N 0 -12/N 0 0 0 0 -4/N 0 0 0 0 0 -10 0 0 0 0 0 0 0 0 0 0 4/N -12/N 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 -16/N*1i 0 0 -8/N*1i -12/N*1i 0 0 0 0 0 0 0 0 0 0 4/N -2/N 0 -10/N -4/N -4/N 0 0 -4/N 0 0 0 0 0 -2 4 0 0 0 0 0 0 0 0 0 0 -24/N -8/N 0 0 2/N*1i 0 0 10*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -2/N 8/N 4/N -6/N -6/N -18/N 0 0 -4/N 0 0 0 0 4 2 0 0 0 0 0 0 0 0 0 0 24/N 8/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N -8/N 0 0 0 0 0 -4/N -12/N -4/N 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -4/N -4/N 0 0 0 0 0 0 0 0 -6/N 0 6/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N -4/N -8/N 0 0 0 0 0 0 0 0 -12/N 0 0 0 0 0 0 -12 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 -8/N*1i 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N -4/N 0 4/N 0 0 -2/N 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 -8/N*1i 0 -6/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 -8*1i 0 0 0 0 0 0 0 0 -12/N*1i -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 4/N -4/N -12/N 0 2/N 6/N 0 4/N 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 -4/N*1i 10/N*1i 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 0 0 0 0 0 2/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 2*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N -8/N 0 12/N 0 0 0 6/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 6*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 -10/N*1i 4/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 -24/N*1i -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N -12/N -4/N 0 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 -6/N -4/N 0 0 0 0 0 -4/N 0 -2 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 -2/N 0 -6/N 0 -6/N 0 -4/N 0 0 -2 0 0 0 0 0 0 -10/N*1i -2/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 2/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 -6*1i 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N 0 6/N 0 8/N 0 -12/N 0 -4/N 0 0 0 0 -2 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 -2*1i 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N -4/N -4/N 0 -6 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 2*1i -2*1i 2*1i 0 -8/N*1i 0 -8/N*1i 0 0 0 0 -8/N*1i 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -14 4/N -20/N 12/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 0 -4*1i -2*1i -2*1i 0 0 0 0 0 0 0 0 -8/N*1i -8/N*1i 0 0 0 -8/N*1i -8/N*1i 0 -8/N*1i 0 0 0 0 0 0 4/N -4 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 -6 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 -4*1i 0 0 0 -4/N*1i 8/N*1i 0 0 0 0 0 8/N*1i 0 0 0 0 0 4/N 0 0 -6 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0 -4/N*1i 0 0 -6/N -4/N 0 -10 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -2/N*1i 6/N*1i 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 -6*1i -2*1i 0 0 0 -8/N*1i -4/N*1i 0 0 -8/N 0 6/N 0 -12];
end
