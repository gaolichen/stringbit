function f = pham6(N)
f=[
    48-60/N 80/N 128/N 144/N 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    48/N 40-60/N 0 0 128/N 96/N 0 0 0 0 0 0 -2/N*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    96/N 0 48-60/N 0 16/N 48/N 192/N 0 0 0 0 2/N*1i 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    48/N 0 0 48-60/N 0 32/N 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 40/N 0 0 32-60/N 0 0 144/N 64/N 0 0 0 0 0 -2/N*1i 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 80/N 32/N 48/N 0 40-60/N 0 48/N 128/N 0 0 0 0 2/N*1i 2/N*1i 0 0 0 2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 32/N 0 0 0 48-60/N 0 16/N 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 32/N 0 0 24-60/N 0 128/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 32/N 24/N 0 0 32-60/N 96/N 0 0 0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 24/N 0 16-60/N 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 -60/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    12*1i 20/N*1i 32/N*1i 36/N*1i 0 0 0 0 0 0 0 56-60/N -4 64/N 0 -4/N 104/N 0 108/N 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -4*1i 0 0 0 0 0 4/N*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -4 56-60/N 16/N 48/N 0 -16/N -4/N -72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0;
    12/N*1i 10*1i 0 0 32/N*1i 24/N*1i 0 0 0 0 0 36/N 12/N 48-60/N -4 0 0 0 0 96/N -4/N 72/N 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -4/N 32/N -4 52-60/N 0 0 0 0 32/N 0 -56/N -4/N 0 0 0 0 0 0 0 0 4/N*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 40-60/N 0 0 0 0 64/N 0 0 96/N 96/N 0 0 0 0 0 -2/N*1i 0 0 -4*1i -2*1i 0 0 0 0 0 0 0 0 0;
    12/N*1i 0 8*1i 0 0 12/N*1i 48/N*1i 0 0 0 0 52/N 4/N 0 0 0 56-60/N 0 0 16/N 0 32/N 0 -4/N 0 0 0 0 0 8/N*1i 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 0 48-60/N 0 0 16/N 0 48/N 0 -48/N 0 0 0 0 0 4/N*1i 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0;
    12/N*1i 0 0 12*1i 0 8/N*1i 0 0 0 0 0 40/N -20/N 0 0 0 0 0 60-60/N 0 0 32/N 0 0 -4/N 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 -4*1i 0 0 0 0 0 0;
    0 10/N*1i 0 0 8*1i 0 0 36/N*1i 16/N*1i 0 0 0 0 28/N 16/N 0 0 0 0 40-60/N 0 0 0 0 0 96/N -4/N 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 32/N 0 0 0 0 32-60/N 0 0 0 0 0 96/N 64/N 0 0 0 0 0 0 -4/N*1i 0 0 0 -4*1i 0 0 0 0 0;
    0 10/N*1i 8/N*1i 0 0 6*1i 0 0 32/N*1i 0 0 0 0 40/N -20/N 0 24/N 0 24/N 0 0 52-60/N 0 0 0 48/N 0 -4/N 0 0 0 0 6/N*1i 0 0 0 0 4/N*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 24/N 0 0 0 0 40-60/N 0 0 0 32/N -32/N 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -6*1i 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48/N -4/N 0 0 0 0 0 0 40-60/N 0 0 16/N 32/N 0 0 0 0 0 2/N*1i 2/N*1i 0 -6/N*1i 0 0 0 -6*1i 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32/N 0 -8/N -4/N 0 0 0 0 0 40-60/N 0 0 32/N 0 0 0 0 0 0 4/N*1i -2/N*1i 0 0 0 0 0 -2*1i 0 0;
    0 0 0 0 8/N*1i 0 0 6*1i 0 32/N*1i 0 0 0 0 0 0 0 0 0 24/N 0 0 0 0 0 36-60/N 0 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 24/N 0 0 0 0 0 24-60/N 0 96/N 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -6*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32/N -4/N -8/N 16/N 24/N 0 0 32-60/N 48/N 0 0 0 0 0 0 0 0 0 4/N*1i -6/N*1i 0 2/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 16/N 0 16-60/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 4*1i 0 8/N*1i 4/N*1i 0 8/N*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 64-60/N 0 0 32/N -4/N 0 0 12/N 32/N 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 64-60/N 0 0 16/N -16/N -8/N 72/N 8/N 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 2*1i 0 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 64-60/N 16/N 0 -4/N 4/N 0 16/N 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 4/N*1i -2/N*1i 4*1i -2*1i 0 0 0 0 16/N*1i 0 8/N*1i 0 0 0 0 0 0 0 24/N 0 24/N 60-60/N 0 0 0 0 0 -4/N 12/N 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 8/N*1i 0 0 8/N*1i 12/N*1i 0 0 0 0 -4/N 8/N 0 0 52-60/N -4 0 0 0 32/N 0 120/N 40/N 0 -10*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 -4/N -8/N 0 -4 48-60/N 0 0 0 32/N 0 -24/N -8/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 -8/N 8/N 0 0 0 56-60/N 0 0 16/N 48/N 0 -16/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 4/N 20/N 0 0 0 0 0 60-60/N 0 0 32/N 32/N 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 4/N*1i -2/N*1i 0 0 0 0 0 6*1i 0 0 8/N*1i 0 0 0 0 0 0 0 24/N 4/N 24/N 0 0 0 0 0 72-60/N 0 0 0 -4/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 4*1i 0 0 0 0 0 16/N*1i 8/N*1i 0 0 0 0 -4/N 16/N 24/N 0 0 0 40-60/N 0 0 0 96/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 4/N 0 0 12/N 16/N 0 0 52-60/N 0 0 32/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 2*1i 0 0 0 4/N*1i 0 0 0 0 0 20/N 0 0 8/N 0 0 0 52-60/N 0 16/N 10/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 -2/N*1i 0 0 0 0 0 0 6*1i 0 0 8/N*1i 0 0 0 0 0 40/N 0 -12/N 0 -8/N 0 0 0 52-60/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 2*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 12/N 0 0 0 36-60/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 2*1i 0 0 0 0 4/N*1i 0 0 4/N*1i 0 60-60/N];
end
