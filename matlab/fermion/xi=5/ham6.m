function f = ham6(N)
f=[
    12 120/N 120/N 192/N 192/N 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 2*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    110/N 10 0 0 0 0 192/N 96/N 144/N 144/N 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 6/N*1i 4*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    22/N 0 10 0 0 0 0 96/N 0 0 144/N 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    96/N 0 0 12 0 0 24/N 0 72/N 0 72/N 192/N 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    48/N 0 0 0 12 0 0 24/N 0 72/N 0 96/N 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    72/N 0 0 0 0 12 0 0 48/N 48/N 48/N 0 0 0 0 0 0 0 0 2/N*1i -2/N*1i -2/N*1i -6/N*1i 0 0 0 0 0 0 6*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 88/N 0 20/N 0 0 8 0 0 0 0 0 216/N 72/N 96/N 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 0;
    0 22/N 110/N 0 20/N 0 0 8 0 0 0 0 0 144/N 0 96/N 0 0 0 0 0 0 0 -2/N*1i 0 0 2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 72/N 0 66/N 0 66/N 0 0 10 0 0 0 72/N 0 96/N 96/N 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 0 2/N*1i 0 0 -2/N*1i 0 0 6*1i 0 0 0 0 0 0 0 0 0 0;
    0 48/N 0 0 88/N 44/N 0 0 0 10 0 0 0 48/N 96/N 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0;
    0 0 120/N 22/N 0 22/N 0 0 0 0 10 0 0 24/N 0 96/N 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -2/N*1i 0 -6/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0;
    0 0 0 48/N 48/N 0 0 0 0 0 0 12 0 0 24/N 24/N 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 66/N 0 20/N 0 0 0 6 0 0 0 192/N 48/N 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 6*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 22/N 88/N 0 20/N 20/N 0 0 6 0 0 0 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 2*1i 0 0 0 0 0;
    0 0 0 0 0 0 48/N 0 44/N 66/N 0 40/N 0 0 8 0 144/N 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 48/N 22/N 0 66/N 20/N 0 0 0 8 0 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -6/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 44/N 0 20/N 0 4 0 120/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 22/N 66/N 20/N 40/N 0 4 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6/N*1i -2/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20/N 20/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -4*1i -8/N*1i 0 -8/N*1i -16/N*1i -12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 0 72/N 0 24/N 0 88/N 40/N 180/N 44/N 0 0 0 0 0 0 0 0 4*1i 0 0 -12/N*1i 0;
    -2*1i -4/N*1i 0 -8/N*1i 0 -6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 4 6 24/N 48/N 24/N -4/N 48/N 32/N 18/N 18/N 0 0 0 0 0 0 0 0 2*1i 0 0 6/N*1i 0;
    -2*1i -4/N*1i 0 -8/N*1i 0 -6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 2 4 0 6 24/N 48/N 24/N 4/N 48/N 24/N 18/N 26/N 0 0 0 0 0 0 0 0 2*1i 0 0 6/N*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 0 0 24/N 0 24/N 8/N 0 0 -16/N 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -4/N*1i -4*1i 0 0 0 0 -16/N*1i 0 -8/N*1i -12/N*1i 0 0 0 0 0 0 0 0 0 62/N 20/N 20/N 0 -2 4 0 0 0 0 0 0 96/N 24/N 120/N 40/N 0 0 0 0 2/N*1i 10*1i 0 0 0;
    -2/N*1i -2*1i 0 0 0 0 -8/N*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 4/N 42/N 42/N 54/N 4 2 0 0 0 0 0 0 96/N 24/N 24/N 8/N 0 0 0 0 0 0 0 0 0;
    -2/N*1i 0 -10*1i 0 0 0 0 -16/N*1i 0 0 -24/N*1i 0 0 0 0 0 0 0 0 22/N 20/N 20/N 0 0 0 2 4 0 0 0 0 0 72/N 0 0 88/N 0 0 0 -4/N*1i 0 4*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -6/N 2/N 66/N 0 0 4 -2 0 0 0 0 0 24/N 0 0 -24/N 4/N 0 0 -2/N*1i 0 -2*1i 0 0;
    -4/N*1i 0 0 -4*1i 0 0 0 0 -8/N*1i 0 0 -16/N*1i 0 0 0 0 0 0 0 48/N 20/N 20/N 0 0 0 0 0 4 0 0 0 24/N 0 72/N 0 24/N 0 0 0 -8/N*1i 0 0 0 0;
    -4/N*1i 0 0 0 -8*1i 0 0 0 0 -12/N*1i 0 -16/N*1i 0 0 0 0 0 0 0 44/N 28/N 20/N -12/N 0 0 0 0 0 4 0 0 0 24/N 0 48/N 0 -4/N 0 0 -4/N*1i 0 0 0 0;
    -2/N*1i 0 0 0 0 -2*1i 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 20/N 0 0 0 0 0 0 0 0 0 0 0 0 0 48/N 0 0 0 0 0 -2/N*1i 0 0 2*1i 0;
    -6/N*1i 0 0 0 0 -6*1i 0 0 -8/N*1i 0 0 0 0 0 0 0 0 0 0 40/N 16/N 24/N -60/N 0 0 0 0 0 0 0 0 0 0 0 48/N 48/N 4/N 0 0 2/N*1i 0 0 6*1i 0;
    0 -4/N*1i 0 0 0 0 -4*1i 0 0 0 0 0 -24/N*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 44/N 36/N 0 0 20/N 0 0 0 0 0 0 0 0 0 216/N 24/N 0 0 0 0 0;
    0 -2/N*1i -10/N*1i 0 0 0 0 -8*1i 0 0 0 0 0 -24/N*1i 0 -16/N*1i 0 0 0 0 0 0 0 22/N 18/N 62/N 14/N 0 20/N 0 0 0 0 0 0 0 0 0 96/N 0 -10/N*1i -2/N*1i 0 0;
    0 -2/N*1i 0 -2/N*1i 0 0 0 0 -2*1i 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0 0 0 20/N 0 0 0 18/N 0 66/N 0 0 0 -2 0 0 0 72/N 0 0 -10/N*1i 0 -2/N*1i 0;
    0 -4/N*1i 0 0 -8/N*1i 0 0 0 0 -6*1i 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 40/N 0 0 0 0 36/N 0 44/N 0 0 0 -2 0 0 0 48/N 0 0 0 0 0;
    0 0 -10/N*1i -2/N*1i 0 0 0 0 0 0 -6*1i 0 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 40/N -20/N 18/N 0 0 22/N 0 0 0 0 -2 0 0 24/N 0 0 -6/N*1i -6/N*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 0 -4/N 0 4/N 0 0 0 0 0 10 0 0 0 0 0 0 6*1i;
    0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 -2*1i 0 0 0 -16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18/N 0 20/N 0 0 0 -6 0 0 0 0 0 0;
    0 0 0 0 0 0 -2/N*1i -8/N*1i 0 0 0 0 0 -6*1i 0 0 0 -24/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 18/N 36/N 0 20/N 20/N 0 0 -6 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i -2*1i -2*1i 0 -8/N*1i -8/N*1i 0 0 -8/N*1i -8/N*1i 0 -8/N*1i 0 0 0 0 0 0 0 0 -4 120/N 24/N 72/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -2*1i 0 0 0 0 0 0 0 -8/N*1i 0 0 -4/N*1i 0 0 0 0 18/N -10 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 6/N*1i 0 0 -4*1i 2*1i 0 0 0 0 0 -8/N*1i 0 0 -8/N*1i 0 0 0 18/N 0 -10 0 -12/N;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -2/N*1i 6/N*1i 0 0 0 0 0 0 -6*1i -2*1i 0 0 -8/N*1i -4/N*1i 0 0 0 0 16/N 0 0 -12 4/N;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 -4/N 4/N -2];
end
