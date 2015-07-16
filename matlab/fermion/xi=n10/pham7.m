function f = pham7(N)
f=[
    140-154/N 240/N 240/N 400/N 400/N 480/N 480/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i -2*1i -2*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    120/N 120-154/N 0 0 0 0 0 400/N 200/N 320/N 320/N 360/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 -4/N*1i -4*1i -2*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    20/N 0 120-154/N 0 0 0 0 0 200/N 0 0 0 320/N 360/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    200/N 0 0 140-154/N 0 0 0 40/N 0 160/N 0 0 160/N 0 480/N 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 2/N*1i -2/N*1i -2/N*1i 0 0 0 0 0 0 -4*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    80/N 0 0 0 140-154/N 0 0 0 40/N 0 160/N 0 0 0 0 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    160/N 0 0 0 0 140-154/N 0 0 0 80/N 0 120/N 0 240/N 0 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    120/N 0 0 0 0 0 140-154/N 0 0 0 80/N 120/N 80/N 0 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 4/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 100/N 0 0 0 0 0 100-154/N 0 0 0 0 0 0 0 0 480/N 160/N 240/N 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -6/N*1i 0 0 0 0 0 0 0 0 0 0 -4*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 20/N 120/N 0 0 0 0 0 100-154/N 0 0 0 0 0 0 0 0 320/N 0 0 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 160/N 0 80/N 0 60/N 0 0 0 120-154/N 0 0 0 0 0 0 120/N 0 240/N 0 120/N 320/N 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 -4/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 80/N 0 0 100/N 0 40/N 0 0 0 120-154/N 0 0 0 0 0 0 80/N 0 240/N 0 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 120/N 0 0 0 60/N 80/N 0 0 0 0 120-154/N 0 0 0 0 0 0 160/N 160/N 80/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 2/N*1i 2/N*1i 6/N*1i 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 -6*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 240/N 20/N 0 0 20/N 0 0 0 0 0 120-154/N 0 0 0 0 40/N 0 0 120/N 0 480/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 0 0 0 0 0 0 6/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 120/N 0 0 20/N 0 0 0 0 0 0 0 120-154/N 0 0 0 0 0 0 80/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 120/N 0 0 80/N 0 0 0 0 0 0 0 140-154/N 0 0 0 40/N 0 0 80/N 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 80/N 200/N 80/N 0 0 0 0 0 0 0 0 0 140-154/N 0 0 0 40/N 40/N 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 2/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 80/N 0 0 0 0 0 0 0 0 80-154/N 0 0 0 0 0 0 480/N 120/N 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 20/N 100/N 0 0 0 0 0 0 0 0 80-154/N 0 0 0 0 0 0 360/N 0 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 120/N 0 60/N 0 60/N 0 0 0 0 0 0 100-154/N 0 0 0 0 240/N 0 240/N 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 -2/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 80/N 0 0 80/N 40/N 0 0 0 0 0 0 0 100-154/N 0 0 0 0 120/N 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 200/N 20/N 0 20/N 80/N 120/N 0 0 0 0 0 0 100-154/N 0 0 0 120/N 0 320/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 2/N*1i 0 6/N*1i 0 0 0 2/N*1i 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 80/N 80/N 0 0 0 40/N 60/N 0 0 0 0 0 120-154/N 0 0 0 120/N 80/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 80/N 0 20/N 0 0 0 0 0 0 0 120-154/N 0 0 0 40/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 60/N 0 0 0 0 0 0 60-154/N 0 0 0 400/N 80/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20/N 80/N 0 0 0 0 0 0 60-154/N 0 0 0 320/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 0 40/N 60/N 0 0 0 0 0 80-154/N 0 400/N 160/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 20/N 0 60/N 0 0 0 0 0 80-154/N 0 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 6/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 0 0 0 40-154/N 0 240/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20/N 60/N 0 0 0 40-154/N 600/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -154/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    4*1i 8/N*1i 0 8/N*1i 20/N*1i 12/N*1i 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 152-154/N -2 -2 0 0 160/N 0 0 0 40/N 0 248/N 0 88/N 0 252/N 84/N 84/N 528/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -2*1i 0 0 0 0 0 0 4/N*1i 12/N*1i 0 0 0 0 0 0 0 0;
    2*1i 4/N*1i 0 8/N*1i 0 6/N*1i 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 152-154/N -2 -2 -2 40/N 120/N 0 0 40/N 4/N -8/N 168/N 32/N 44/N 126/N 84/N 42/N -24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i -2*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0;
    2*1i 4/N*1i 0 8/N*1i 0 6/N*1i 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 152-154/N -2 -2 40/N 0 120/N 0 40/N -4/N -8/N 168/N 40/N 44/N 126/N 84/N 34/N -24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i -2*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0;
    2*1i 4/N*1i 0 8/N*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 152-154/N -4 0 80/N 80/N 0 40/N 0 176/N -8/N 80/N -88/N -12/N 36/N 92/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 -4/N*1i -12/N*1i 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -4 148-154/N 0 40/N 40/N 240/N 0 80/N -8/N 72/N 0 80/N -12/N -48/N -92/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0;
    4/N*1i 4*1i 0 0 0 0 0 16/N*1i 0 8/N*1i 16/N*1i 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 84/N 24/N 24/N -4/N 0 132-154/N -2 -2 0 0 0 0 0 0 0 0 0 0 0 240/N 0 40/N 0 168/N 88/N 396/N 84/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 -4*1i 0 0 0 0 0 0 0 0 0 12/N*1i 0 0 0 0 0;
    2/N*1i 2*1i 0 0 0 0 0 8/N*1i 0 8/N*1i 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N 64/N -2/N 42/N 22/N -2 132-154/N -4 -6 0 0 0 0 0 0 0 0 0 0 80/N 160/N 40/N 4/N 80/N 32/N -18/N 46/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2*1i 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0;
    2/N*1i 2*1i 0 0 0 0 0 8/N*1i 0 8/N*1i 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N -2/N 64/N 42/N 22/N -2 -4 132-154/N -6 0 0 0 0 0 0 0 0 0 0 80/N 160/N 40/N -4/N 80/N 40/N -18/N 38/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2*1i 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N -2/N 0 44/N 0 -2 -2 132-154/N 0 0 0 0 0 0 0 0 0 0 0 80/N 0 40/N -8/N 0 0 -48/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    2/N*1i 0 12*1i 0 0 0 0 0 20/N*1i 0 0 0 32/N*1i 36/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20/N 22/N 22/N 20/N 0 0 0 0 0 128-154/N -4 0 0 0 0 0 0 0 0 0 0 160/N 0 0 0 0 0 248/N 0 252/N 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 2/N*1i 0 -4*1i 0 -4*1i 0 0 0 0 0 0 0 0 0 4/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 4/N -4/N 4/N 40/N 0 0 0 0 -4 128-154/N 0 0 0 0 0 0 0 0 0 0 40/N 120/N 0 0 0 0 -16/N -4/N -144/N 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0;
    4/N*1i 0 0 4*1i 0 0 0 0 0 8/N*1i 0 0 0 0 16/N*1i 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 124/N -2/N -2/N 88/N 0 0 0 0 0 0 0 152-154/N -4 0 0 0 0 0 0 40/N 0 0 0 80/N 0 0 0 40/N 0 0 528/N 88/N 0 0 0 0 0 0 0 0 6/N*1i 0 2/N*1i 0 0 0 0 -10*1i 0 0 0 0 0 0 0 0 0 0 0;
    2/N*1i 0 0 2*1i 0 0 0 0 0 4/N*1i 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 82/N 82/N 0 44/N 0 0 0 0 0 0 -4 148-154/N 0 0 0 0 0 0 0 40/N 0 0 80/N 0 0 0 40/N 0 0 -48/N -8/N 0 0 0 0 0 0 0 0 0 6/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    4/N*1i 0 0 0 10*1i 0 0 0 0 0 16/N*1i 0 0 0 0 24/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 84/N 34/N 42/N 84/N 4/N 0 0 0 0 0 0 0 0 148-154/N -4 0 0 0 0 0 0 40/N 0 0 120/N 0 0 0 4/N 0 0 168/N 0 0 0 0 0 0 0 0 0 4/N*1i 4/N*1i 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 42/N 42/N -84/N 84/N 0 0 0 0 0 0 0 0 -4 152-154/N 0 0 0 0 0 0 0 40/N 0 40/N 0 0 0 0 0 0 -104/N 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0;
    4/N*1i 0 0 0 0 4*1i 0 0 0 0 0 8/N*1i 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 44/N 44/N 0 0 0 0 0 0 0 0 0 0 0 0 148-154/N 0 0 0 0 0 0 0 80/N 0 120/N 0 0 0 40/N 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 4/N*1i 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0;
    6/N*1i 0 0 0 0 0 8*1i 0 0 0 0 12/N*1i 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 84/N 84/N 84/N 36/N -36/N 0 0 0 0 0 0 0 0 0 0 0 148-154/N 0 0 0 0 0 0 0 80/N 0 80/N 80/N 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i -2/N*1i 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0;
    8/N*1i 0 0 0 0 6*1i 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 88/N 48/N 40/N 88/N -88/N 0 0 0 0 0 0 0 0 0 0 0 0 152-154/N 0 0 0 0 0 0 0 0 120/N 0 -4/N 120/N 0 160/N 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0;
    2/N*1i 0 0 0 0 0 2*1i 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 44/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 152-154/N 0 0 0 0 0 0 120/N 0 0 0 0 160/N 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0;
    0 4/N*1i 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 24/N*1i 0 8/N*1i 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 64/N 22/N 22/N 0 0 0 0 0 0 0 0 0 0 0 112-154/N -4 0 0 0 0 0 0 0 0 0 0 0 240/N 40/N 264/N 88/N 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 -10*1i 0 0 0 0 0 0 0;
    0 2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 12/N*1i 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 42/N 42/N 72/N 0 0 0 0 0 0 0 0 0 0 -4 108-154/N 0 0 0 0 0 0 0 0 0 0 0 240/N 40/N -24/N -8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 2/N*1i 12/N*1i 0 0 0 0 0 10*1i 0 0 0 0 0 0 0 0 32/N*1i 0 0 24/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20/N 22/N 22/N 0 84/N 24/N 0 0 0 0 0 0 0 0 0 0 108-154/N -4 0 0 0 0 0 0 0 0 0 0 240/N 0 0 168/N 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 6/N -2/N 60/N -4/N 68/N 0 0 0 0 0 0 0 0 0 0 -4 112-154/N 0 0 0 0 0 0 0 0 0 0 80/N 0 0 -104/N -4/N 0 0 0 0 0 2/N*1i 0 0 4/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0;
    0 4/N*1i 0 4/N*1i 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 16/N*1i 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 44/N 44/N 0 0 0 40/N 48/N 0 0 60/N 0 0 0 0 0 0 0 128-154/N 0 0 0 0 0 0 0 0 120/N 0 240/N 0 40/N 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0;
    0 4/N*1i 0 0 10/N*1i 0 0 0 0 0 8*1i 0 0 0 0 0 0 0 0 24/N*1i 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 84/N 36/N 44/N 12/N 0 0 0 0 64/N 28/N 0 40/N 0 0 0 0 0 0 0 128-154/N 0 0 0 0 0 0 0 0 80/N 0 160/N 0 4/N 0 0 0 0 0 4/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0;
    0 2/N*1i 0 0 0 2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 44/N 0 0 0 0 0 0 0 0 0 24/N 0 0 80/N 0 0 0 0 0 0 132-154/N 0 0 0 0 0 0 0 0 160/N 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0;
    0 6/N*1i 0 0 0 0 8/N*1i 0 0 0 0 6*1i 0 0 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 88/N 48/N 40/N -132/N 0 0 0 0 0 0 0 48/N 60/N 0 0 0 0 0 0 0 0 132-154/N 0 0 0 0 0 0 0 0 160/N 80/N -4/N 0 0 0 0 0 -2/N*1i 0 0 0 0 0 -2/N*1i 0 0 0 -6*1i 0 0 0 0 0;
    0 0 12/N*1i 2/N*1i 0 0 0 0 0 0 0 0 8*1i 0 0 0 0 0 0 0 12/N*1i 0 48/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 124/N 4/N 20/N 24/N 0 0 0 20/N 0 0 0 0 0 0 0 0 0 0 128-154/N 0 0 0 0 0 40/N 0 0 80/N 0 0 0 0 0 0 0 8/N*1i 0 4/N*1i 10/N*1i 0 0 6/N*1i 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 4/N 0 0 0 -4/N 0 0 0 0 0 0 0 0 0 0 120-154/N 0 0 0 0 0 0 0 0 120/N 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0;
    0 0 12/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 12*1i 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 88/N -44/N 0 0 0 0 24/N 0 20/N 0 0 0 0 0 0 0 0 0 0 0 132-154/N 0 0 0 0 0 0 80/N 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 2/N*1i 0 0 0 0 0 0 -4*1i 0 0;
    0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 44/N 0 0 0 0 0 0 80/N 0 0 0 0 0 0 0 0 0 0 0 152-154/N 0 0 0 40/N 0 0 0 0 0 0 0 0 0 0 0 0 10/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0;
    0 0 0 4/N*1i 10/N*1i 0 0 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 88/N 0 88/N -44/N 0 0 80/N 0 0 0 0 0 0 0 0 0 0 0 0 0 152-154/N 0 0 0 40/N 40/N 0 0 0 0 0 0 0 0 0 0 0 6/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 32/N*1i 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 48/N 0 0 0 0 0 0 0 0 0 0 0 88-154/N 0 0 0 0 0 480/N 40/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 2/N*1i 10/N*1i 0 0 0 0 0 0 0 0 8*1i 0 0 0 0 0 0 36/N*1i 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20/N 24/N 64/N 28/N 0 0 0 0 0 0 0 0 0 0 88-154/N 0 0 0 0 0 240/N 0 0 0 0 0 0 0 0 0 0 0 10/N*1i 2/N*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 44/N 0 0 0 24/N 0 60/N 0 0 0 0 0 0 0 0 112-154/N 0 0 0 240/N 0 0 0 0 0 0 0 0 0 0 0 0 10/N*1i 0 2/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 24/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 88/N 0 0 0 0 48/N 0 40/N 0 0 0 0 0 0 0 0 112-154/N 0 0 0 120/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 10/N*1i 2/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 32/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 88/N -44/N 24/N 0 0 20/N 48/N 0 60/N 0 0 0 0 0 0 112-154/N 0 0 120/N 0 0 0 0 0 0 0 0 0 0 0 0 6/N*1i 6/N*1i 0 0 4/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 4/N 0 -4/N 0 60/N 0 0 0 0 0 0 0 0 100-154/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 -6*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 2*1i 0 0 0 20/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 24/N 0 0 0 0 0 72-154/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 8/N*1i 0 0 0 0 0 0 6*1i 0 0 0 32/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 24/N 48/N 0 0 0 0 0 72-154/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 4*1i 0 8/N*1i 4/N*1i 4/N*1i 0 0 0 8/N*1i 0 8/N*1i -4/N*1i 12/N*1i 4/N*1i 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 160-154/N -4 0 80/N 40/N 0 0 440/N 88/N 88/N 132/N 0 0 0 0 0 0 0 -14*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 2*1i 0 2*1i 0 4/N*1i 4/N*1i 12/N*1i 0 0 0 8/N*1i 0 8/N*1i 0 4/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4 160-154/N -4 80/N 40/N 4/N 0 -40/N -8/N -8/N -12/N 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 2*1i 2*1i 0 0 4/N*1i 4/N*1i 4/N*1i 0 0 0 8/N*1i 8/N*1i 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4 156-154/N 80/N 0 -4/N 80/N 0 0 80/N 120/N 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 2/N*1i 0 0 4*1i 2*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 16/N*1i 16/N*1i 0 0 8/N*1i 8/N*1i 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 48/N 48/N 136-154/N 0 0 0 0 0 0 0 400/N 40/N 120/N 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 -2/N*1i 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 4/N*1i 0 0 0 0 8/N*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 20/N 20/N 4/N 0 136-154/N 0 0 0 0 0 0 0 80/N 0 0 12/N 80/N 0 14/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N 0 0 136-154/N 0 0 0 0 0 0 0 0 -8/N 144/N 8/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 20/N 0 0 0 136-154/N 0 0 0 0 0 40/N 0 4/N 0 40/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 44/N 0 0 0 0 0 0 160-154/N 0 0 0 40/N 0 0 0 0 0 0 14/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 4/N*1i -4/N*1i 0 0 0 0 0 0 0 0 4*1i -2*1i 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 88/N 0 0 0 0 0 0 0 160-154/N 0 0 0 40/N 0 4/N 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 4/N*1i 4/N*1i 4/N*1i -4/N*1i 0 0 0 0 0 0 0 0 0 0 6*1i 0 4*1i 0 0 0 0 0 8/N*1i 0 0 8/N*1i 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 88/N 8/N 88/N 0 0 0 0 0 0 160-154/N 0 0 0 120/N -4/N 0 80/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 2/N*1i 2/N*1i -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 8*1i 0 0 0 0 0 4/N*1i 12/N*1i 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 44/N 4/N 44/N 0 0 0 0 0 0 0 160-154/N 0 0 80/N 0 -4/N 0 0 14/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 12/N*1i 0 0 4/N*1i 0 0 0 0 0 0 0 24/N 0 0 0 0 0 0 0 120-154/N 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i -6/N*1i 4/N*1i -2/N*1i 0 0 0 0 0 0 0 0 0 0 4*1i -2*1i 0 0 0 0 0 0 0 0 0 0 16/N*1i 0 0 8/N*1i 0 0 0 0 0 0 24/N 48/N 0 48/N 0 0 0 0 0 120-154/N 0 0 0 0 12/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 2/N*1i -6/N*1i 0 0 0 0 0 0 0 0 2/N*1i 8/N*1i 0 0 0 0 0 0 6*1i 2*1i 0 0 0 0 0 0 0 16/N*1i 8/N*1i 0 0 0 0 0 0 0 48/N 0 0 0 0 0 24/N 48/N 0 0 144-154/N 0 0 0 -4/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 -8/N 8/N 0 4/N -4/N 0 0 0 0 128-154/N 0 0 120/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 4/N 44/N 0 0 0 0 -4/N 0 0 0 0 132-154/N 0 80/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i -2/N*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 48/N 4/N 48/N 0 0 24/N 0 0 0 0 0 0 144-154/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N 24/N 40/N 0 112-154/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 4/N*1i 0 0 0 0 4/N*1i 4/N*1i 0 0 0 0 0 0 0 0 168-154/N];
end
