function f = pham6(N)
f=[
    120-132/N 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -120/N 100-132/N 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 120-132/N 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 120-132/N 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -100/N -40/N 0 80-132/N 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 -80/N -120/N 0 100-132/N 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 0 2/N*1i 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 120-132/N 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 -80/N -40/N 0 60-132/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 -60/N -120/N 0 80-132/N 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 -60/N -80/N 40-132/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 -40/N -132/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    12*1i 20/N*1i 32/N*1i 36/N*1i 0 0 0 0 0 0 0 128-132/N -4 0 0 -4/N 8/N 0 12/N 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -4*1i 0 0 0 0 0 4/N*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -4 128-132/N 0 0 0 -16/N -4/N -24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0 0 0 0 0;
    12/N*1i 10*1i 0 0 32/N*1i 24/N*1i 0 0 0 0 0 -76/N -16/N 108-132/N -4 0 0 0 0 0 -4/N 8/N 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -4/N -52/N -4 112-132/N 0 0 0 0 0 0 -24/N -4/N 0 0 0 0 0 0 0 0 4/N*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 100-132/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 -4*1i -2*1i 0 0 0 0 0 0 0 0 0;
    12/N*1i 0 8*1i 0 0 12/N*1i 48/N*1i 0 0 0 0 4/N 4/N 0 0 0 128-132/N 0 0 0 0 0 0 -4/N 0 0 0 0 0 8/N*1i 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 0 120-132/N 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0;
    12/N*1i 0 0 12*1i 0 8/N*1i 0 0 0 0 0 8/N -4/N 0 0 0 0 0 132-132/N 0 0 0 0 0 -4/N 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 -4*1i 0 0 0 0 0 0;
    0 10/N*1i 0 0 8*1i 0 0 36/N*1i 16/N*1i 0 0 0 0 -56/N -12/N 0 -40/N 0 0 88-132/N 0 0 0 0 0 0 -4/N 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 -80/N 0 -40/N 0 0 80-132/N 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 -4*1i 0 0 0 0 0;
    0 10/N*1i 8/N*1i 0 0 6*1i 0 0 32/N*1i 0 0 0 0 8/N -4/N 0 -32/N 0 -60/N 0 0 112-132/N 0 0 0 0 0 -4/N 0 0 0 0 6/N*1i 0 0 0 0 4/N*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 -60/N 0 0 0 0 100-132/N 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -6*1i 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 0 0 0 0 100-132/N 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 -6/N*1i 0 0 0 -6*1i 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20/N -4/N 0 0 0 0 0 100-132/N 0 0 0 0 0 0 0 0 0 4/N*1i -2/N*1i 0 0 0 0 0 -2*1i 0 0;
    0 0 0 0 8/N*1i 0 0 6*1i 0 32/N*1i 0 0 0 0 0 0 0 0 0 -32/N 0 -40/N 0 0 0 72-132/N 0 0 -4/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N -60/N 0 -40/N -40/N 0 0 60-132/N 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -6*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 20/N -40/N -60/N 0 0 80-132/N 0 0 0 0 0 0 0 0 0 0 4/N*1i -6/N*1i 0 2/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N -40/N -40/N 40-132/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 4*1i 0 8/N*1i 4/N*1i 0 8/N*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 136-132/N 0 0 0 -4/N 0 0 12/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 136-132/N 0 0 0 0 -8/N 24/N 8/N 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 2*1i 0 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 136-132/N 0 0 -4/N 4/N 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 4/N*1i -2/N*1i 4*1i -2*1i 0 0 0 0 16/N*1i 0 8/N*1i 0 0 0 0 0 0 0 -32/N 0 -32/N 120-132/N 0 0 0 0 0 -4/N 12/N 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 8/N*1i 0 0 8/N*1i 12/N*1i 0 0 0 0 -4/N -20/N 0 0 112-132/N -4 0 0 0 0 0 24/N 8/N 0 -10*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 24/N -8/N 0 -4 108-132/N 0 0 0 0 0 -24/N -8/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 -8/N 8/N 0 0 0 128-132/N 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 4/N 4/N 0 0 0 0 0 132-132/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 4/N*1i -2/N*1i 0 0 0 0 0 6*1i 0 0 8/N*1i 0 0 0 0 0 0 0 8/N 4/N 8/N 0 0 0 0 0 144-132/N 0 0 0 -4/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 4*1i 0 0 0 0 0 16/N*1i 8/N*1i 0 0 0 0 -4/N -40/N -32/N -40/N 0 0 88-132/N 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 4/N 0 0 -16/N -40/N 0 0 112-132/N 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 2*1i 0 0 0 4/N*1i 0 0 0 0 0 4/N 0 0 -20/N 0 0 0 112-132/N 0 0 10/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 -2/N*1i 0 0 0 0 0 0 6*1i 0 0 8/N*1i 0 0 0 0 0 8/N 0 16/N 0 -8/N 0 0 0 112-132/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 2*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 -16/N -40/N -40/N 0 72-132/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 2*1i 0 0 0 0 4/N*1i 0 0 4/N*1i 0 120-132/N];
end
