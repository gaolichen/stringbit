function f = ham7(N)
f=[
    14 48/N 80/N 96/N 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    42/N 12 0 0 80/N 64/N 72/N 0 0 0 0 0 0 0 0 0 2/N*1i 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    56/N 0 14 0 8/N 32/N 0 96/N 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    56/N 0 0 14 0 16/N 48/N 32/N 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 0 0 0 0 2*1i 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 36/N 4/N 0 10 0 0 0 96/N 48/N 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 48/N 30/N 18/N 0 12 0 0 24/N 48/N 96/N 0 0 0 0 0 0 0 -2/N*1i 0 0 0 2/N*1i 0 0 0 -2/N*1i 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 24/N 0 24/N 0 0 12 0 0 32/N 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 40/N 16/N 0 0 0 14 0 8/N 48/N 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 30/N 4/N 0 0 8 0 0 96/N 32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 40/N 24/N 36/N 8/N 0 10 0 48/N 96/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 16/N 0 18/N 0 0 12 0 24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 24/N 4/N 0 6 0 80/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 16/N 18/N 12/N 0 8 80/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 18/N 8/N 4 168/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14*1i -24/N*1i -40/N*1i -48/N*1i 0 0 0 0 0 0 0 0 0 0 0 6 4 0 40/N 0 4/N 56/N 0 0 60/N 0 48/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 6 4 8/N 32/N 0 16/N 40/N 4/N 12/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 2*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 2 0 16/N 0 8/N -8/N 0 0 4/N -16/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i -2*1i 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i -12*1i 0 0 -40/N*1i -32/N*1i -36/N*1i 0 0 0 0 0 0 0 0 26/N 2/N 4/N 4 4 0 0 0 0 0 0 0 64/N 0 4/N 40/N 0 36/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 2/N*1i 4*1i 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 20/N 4/N 4 4 0 0 0 0 0 0 0 16/N 48/N 0 16/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 12 0 0 0 0 0 0 0 0 40/N 0 0 0 64/N 64/N 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 -2/N*1i 0 0 0 0 4*1i 2*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i 0 -10*1i 0 0 -16/N*1i 0 -48/N*1i 0 0 0 0 0 0 0 28/N 0 0 0 0 0 6 4 0 0 0 0 8/N 0 0 24/N 0 0 4/N 0 0 48/N 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i -2/N*1i -4/N*1i -2/N*1i 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 16/N -12/N 0 0 0 4 2 0 0 0 0 0 8/N 0 8/N 0 0 0 0 0 16/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i -2/N*1i -4/N*1i -2/N*1i 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 14 0 0 0 0 0 8/N 0 32/N 0 0 -32/N 0 0 48/N 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 4*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i 0 0 -8*1i 0 0 -24/N*1i -16/N*1i 0 0 0 0 0 0 0 20/N -4/N -4/N 0 0 0 0 0 0 6 0 0 0 0 0 16/N 0 16/N 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -2/N*1i -2/N*1i -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 14 0 0 0 0 0 16/N 0 16/N 0 -24/N 0 -32/N 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -14/N*1i 0 0 -6*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 0 8/N -4/N 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 24/N 0 4/N 0 32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i -2/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -12/N*1i 0 0 -10*1i 0 0 0 -48/N*1i -24/N*1i 0 0 0 0 0 0 0 0 20/N 2/N 0 4/N 0 0 0 0 0 2 4 0 0 0 0 0 0 0 0 0 72/N 4/N 24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 10/N 0 0 4/N 0 0 0 0 4 -2 0 0 0 0 0 0 0 0 0 24/N 0 8/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 30/N 0 0 4/N 0 0 0 0 0 10 0 0 0 0 0 0 0 0 0 64/N 0 0 48/N 48/N 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 2/N*1i 2/N*1i 6/N*1i 0 0 0 0 0 0 0 4*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -12/N*1i -10/N*1i 0 0 -8*1i 0 0 0 -24/N*1i -48/N*1i 0 0 0 0 0 0 0 20/N -4/N 0 14/N -2/N 0 18/N 0 0 0 0 0 4 0 0 0 0 0 0 0 24/N 0 32/N 0 4/N 0 0 0 0 0 0 0 0 0 0 0 -8/N*1i 0 -4/N*1i 0 0 0 0 -2/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 24/N 0 12/N 0 0 0 0 0 12 0 0 0 0 0 0 0 16/N 0 48/N 0 -24/N 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -12/N*1i 0 -8/N*1i 0 0 -12*1i 0 0 -16/N*1i 0 0 0 0 0 0 0 0 8/N -4/N 0 0 0 0 4/N 0 24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 32/N 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32/N 4/N 0 0 0 6/N 0 0 0 0 0 0 0 12 0 0 0 0 0 8/N 0 0 24/N 0 64/N 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16/N 0 0 -6/N 0 0 4/N 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 24/N 32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 24/N 0 0 0 4/N -6/N 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 16/N 16/N 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i -2/N*1i -2/N*1i -6/N*1i 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 6*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 -10/N*1i 0 0 0 0 -6*1i 0 0 -24/N*1i 0 0 0 0 0 0 0 0 0 0 8/N -4/N 0 0 0 16/N 0 0 0 0 0 0 0 0 0 2 0 0 0 8/N 0 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 24/N 0 -16/N 0 0 0 0 0 0 0 0 0 0 0 14 0 0 0 8/N 8/N 0 -32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 -10/N*1i 0 0 0 -8*1i 0 0 -48/N*1i -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 14/N -2/N 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 64/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 24/N 0 4/N 0 4/N 0 0 0 0 0 8 0 0 0 0 0 0 72/N 32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 -10/N*1i -8/N*1i 0 0 0 -6*1i 0 0 -48/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N -4/N 0 4/N 0 18/N 0 0 0 8/N 0 0 0 -2 0 0 0 0 48/N 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 18/N 0 0 0 0 0 4/N 0 0 0 10 0 0 0 0 24/N -16/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 24/N 4/N 0 0 18/N 0 18/N 0 4/N 0 0 0 0 10 0 0 0 24/N 32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 6/N*1i 0 2/N*1i 0 0 -2/N*1i 0 0 0 6*1i 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16/N 0 -6/N 4/N 0 24/N 12/N 0 0 0 0 0 0 0 10 0 0 0 32/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16/N 16/N 0 4/N -6/N 0 0 0 0 0 0 12 0 0 8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i -2/N*1i 0 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -8/N*1i 0 0 -6*1i 0 -40/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 4/N 0 0 0 0 -6 0 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 18/N 0 4/N 4/N 0 0 0 6 0 64/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 6*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16/N 4/N -6/N 12/N 18/N 8/N 0 0 8 48/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 6/N*1i 0 -2/N*1i 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 12/N 4/N 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -2*1i -8/N*1i -4/N*1i 0 -8/N*1i 0 0 -12/N*1i 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 2 0 2 0 24/N 0 0 4/N 0 0 0 24/N 0 0 0 12/N 12/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 2*1i 0 -4/N*1i 0 0 -8/N*1i 4/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 2 -2 -2 8/N 16/N 0 8/N -8/N 0 0 0 16/N -4/N 2/N 4/N -6/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 6/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 -8/N*1i 0 0 -16/N*1i 0 0 -12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 2 4 8/N 0 32/N 0 4/N 4/N 0 8/N -4/N 0 -4/N 20/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 2*1i 0 -4/N*1i 0 0 -8/N*1i -4/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 2 -2 -2 8/N -16/N 0 -8/N 0 8/N 0 0 -8/N 4/N -10/N 12/N 6/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 -6/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 -4/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 4 -2 -6 0 0 16/N 0 8/N -8/N 12/N 16/N 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 2/N*1i -4*1i 0 0 0 0 0 0 0 0 -16/N*1i -8/N*1i 0 -8/N*1i 0 -12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 14/N 0 2/N 0 0 -4 0 0 0 0 0 0 0 0 0 0 0 0 32/N 4/N 0 0 -12/N 16/N 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N 0 -8/N 0 0 -4 0 0 0 0 0 0 0 0 0 0 0 0 8/N -8/N 8/N 0 -8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 -8/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N 8/N -2/N 2/N 0 0 -4 0 0 0 0 0 0 0 0 0 0 16/N 0 4/N -4/N 0 8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 -8/N*1i 0 0 0 -8/N*1i -16/N*1i -12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 6/N 0 -6/N 0 0 0 0 0 2 2 0 0 0 0 0 0 0 0 24/N 0 0 0 0 24/N 8/N 36/N 12/N 0 0 0 0 0 0 0 4*1i 0 0 0 -12/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 -8/N*1i 0 -6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N -6/N 2/N 2/N 8/N 0 0 0 2 0 4 6 0 0 0 0 0 0 0 8/N 16/N 0 0 0 16/N 16/N 18/N 2/N 0 0 0 0 0 0 -2/N*1i 2*1i 0 0 0 6/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 -8/N*1i 0 -6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N -2/N 6/N 6/N -8/N 0 0 0 2 4 0 6 0 0 0 0 0 0 0 8/N 16/N 0 0 0 16/N 8/N 18/N 10/N 0 0 0 0 0 0 2/N*1i 2*1i 0 0 0 6/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N 0 2/N 4/N 0 0 0 0 2 2 0 0 0 0 0 0 0 0 0 8/N 0 0 0 8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 -4*1i 2*1i 0 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 8/N -4/N 4/N -4/N 0 0 0 0 0 0 0 0 -6 0 0 0 0 0 8/N 0 0 0 0 0 4/N 0 0 0 -12/N 0 0 0 0 0 -10/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 14/N -4/N -6/N 0 0 0 0 0 0 0 0 0 2 4 0 0 0 0 8/N 0 16/N 0 0 0 -8/N 0 0 24/N 0 0 0 0 0 0 0 10*1i 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -6/N 0 6/N -4/N 0 0 0 0 0 0 0 0 4 6 0 0 0 0 0 8/N 16/N 0 0 0 -8/N 0 0 24/N 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 8/N -4/N -16/N 0 0 0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 16/N 24/N 0 16/N 0 0 -8/N 0 0 0 0 0 0 0 0 0 6*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -10/N*1i 2/N*1i 2/N*1i 0 0 0 0 0 0 -6*1i 0 -8*1i 0 0 0 -8/N*1i 0 -12/N*1i 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 4/N -8/N 12/N 0 -12/N 0 0 0 0 0 0 0 0 0 0 0 -6 0 0 0 0 0 0 32/N 0 4/N 0 4/N 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N 0 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 -24/N 0 0 0 -24/N 0 -32/N 0 0 0 0 0 0 0 0 4*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 0 0 0 -4*1i 2*1i 0 0 0 0 0 0 0 0 0 -24/N*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 4/N 0 0 0 0 4/N 0 0 0 0 0 -10 0 0 0 0 0 0 0 0 0 0 4/N -12/N 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 -16/N*1i 0 0 -8/N*1i -12/N*1i 0 0 0 0 0 0 0 0 0 0 4/N 6/N 0 14/N 4/N 4/N 0 0 4/N 0 0 0 0 0 -2 4 0 0 0 0 0 0 0 0 32/N 0 24/N 8/N 0 0 2/N*1i 0 0 10*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -10/N 8/N 4/N 10/N 10/N 6/N 0 0 4/N 0 0 0 0 4 2 0 0 0 0 0 0 0 0 32/N 0 24/N 8/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N -8/N 0 0 0 0 0 12/N 4/N 12/N 0 0 0 0 0 4 0 0 0 0 0 0 0 16/N 48/N 0 -8/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 4/N 0 0 0 0 0 0 0 0 2/N 0 -18/N 0 0 0 0 0 0 0 0 0 0 0 0 32/N 16/N 0 0 0 0 0 -2/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 2/N*1i 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 -12 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 -8/N*1i 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 16/N 4/N 4/N 0 4/N 0 0 6/N 0 0 0 0 0 0 0 0 4 0 0 0 0 8/N 0 24/N 0 0 0 -8/N*1i 0 -6/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 -8*1i 0 0 0 0 0 0 0 0 -12/N*1i -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 12/N 12/N 4/N -12/N 0 -6/N -2/N 0 4/N 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 16/N 0 0 -4/N*1i 10/N*1i 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 0 0 0 0 0 0 -6/N 0 0 0 0 0 0 0 0 0 0 0 0 0 16/N 0 0 0 -2/N*1i 0 2/N*1i 0 2*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N 0 8/N -12/N 0 0 0 -2/N 4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16/N 0 0 2/N*1i 0 0 0 6*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 4/N 0 0 0 -16/N 0 0 0 0 0 0 0 0 0 0 2 0 8/N 8/N 0 0 0 0 -10/N*1i 4/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 -24/N*1i -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 12/N 4/N 4/N 0 0 4/N 0 0 0 0 0 0 0 0 72/N 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 2/N 12/N 0 0 0 0 0 4/N 0 -2 0 0 24/N 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 0 6/N 0 2/N 0 18/N 0 4/N 0 0 -2 0 24/N 0 0 0 0 -10/N*1i -2/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 2/N*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 -6*1i 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8/N 0 -2/N 0 8/N 0 4/N 0 12/N 0 0 0 0 -2 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 -2*1i 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N 4/N 4/N 0 -6 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 2*1i -2*1i 2*1i 0 -8/N*1i 0 -8/N*1i 0 0 0 0 -8/N*1i 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -14 4/N -20/N 12/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 0 -4*1i -2*1i -2*1i 0 0 0 0 0 0 0 0 -8/N*1i -8/N*1i 0 0 0 -8/N*1i -8/N*1i 0 -8/N*1i 0 0 0 0 0 0 4/N -4 0 0 40/N 24/N;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 -6 0 8/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 -4*1i 0 0 0 -4/N*1i 8/N*1i 0 0 0 0 0 8/N*1i 0 0 0 0 0 4/N 0 0 -6 0 -8/N;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0 -4/N*1i 0 0 2/N 4/N 0 -10 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -2/N*1i 6/N*1i 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 -6*1i -2*1i 0 0 0 -8/N*1i -4/N*1i 0 0 0 0 -2/N 0 -12];
end
