function f = pham7(N)
f=[
    -14 216/N 216/N 360/N 360/N 432/N 432/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i -2*1i -2*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    228/N -12 0 0 0 0 0 360/N 180/N 288/N 288/N 324/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i 0 -4/N*1i -4*1i -2*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    38/N 0 -12 0 0 0 0 0 180/N 0 0 0 288/N 324/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    180/N 0 0 -14 0 0 0 36/N 0 144/N 0 0 144/N 0 432/N 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 2/N*1i -2/N*1i -2/N*1i 0 0 0 0 0 0 -4*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    72/N 0 0 0 -14 0 0 0 36/N 0 144/N 0 0 0 0 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    144/N 0 0 0 0 -14 0 0 0 72/N 0 108/N 0 216/N 0 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    108/N 0 0 0 0 0 -14 0 0 0 72/N 108/N 72/N 0 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 4/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 190/N 0 40/N 0 0 0 -10 0 0 0 0 0 0 0 0 432/N 144/N 216/N 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2/N*1i -6/N*1i 0 0 0 0 0 0 0 0 0 0 -4*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 38/N 228/N 0 40/N 0 0 0 -10 0 0 0 0 0 0 0 0 288/N 0 0 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 144/N 0 152/N 0 114/N 0 0 0 -12 0 0 0 0 0 0 108/N 0 216/N 0 108/N 288/N 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 -4/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 72/N 0 0 190/N 0 76/N 0 0 0 -12 0 0 0 0 0 0 72/N 0 216/N 0 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 108/N 0 0 0 114/N 152/N 0 0 0 0 -12 0 0 0 0 0 0 144/N 144/N 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 2/N*1i 2/N*1i 6/N*1i 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 -6*1i -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 216/N 38/N 0 0 38/N 0 0 0 0 0 -12 0 0 0 0 36/N 0 0 108/N 0 432/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 0 0 0 0 0 0 6/N*1i 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 108/N 0 0 38/N 0 0 0 0 0 0 0 -12 0 0 0 0 0 0 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 108/N 0 0 72/N 0 0 0 0 0 0 0 -14 0 0 0 36/N 0 0 72/N 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 72/N 180/N 72/N 0 0 0 0 0 0 0 0 0 -14 0 0 0 36/N 36/N 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 2/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 152/N 0 40/N 0 0 0 0 0 0 -8 0 0 0 0 0 0 432/N 108/N 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 38/N 190/N 0 40/N 0 40/N 0 0 0 0 -8 0 0 0 0 0 0 324/N 0 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 -2/N*1i 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 108/N 0 114/N 0 114/N 0 0 80/N 0 0 0 -10 0 0 0 0 216/N 0 216/N 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 -2/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 72/N 0 0 152/N 76/N 0 0 0 40/N 0 0 0 -10 0 0 0 0 108/N 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 180/N 38/N 0 38/N 152/N 228/N 0 40/N 0 0 0 0 -10 0 0 0 108/N 0 288/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 2/N*1i 0 6/N*1i 0 0 0 2/N*1i 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 72/N 72/N 0 0 0 76/N 114/N 0 0 0 0 0 -12 0 0 0 108/N 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 72/N 0 38/N 0 0 0 0 0 0 0 -12 0 0 0 36/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 114/N 0 40/N 0 0 0 0 -6 0 0 0 360/N 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 38/N 152/N 0 40/N 40/N 0 0 0 -6 0 0 0 288/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 72/N 0 76/N 114/N 0 80/N 0 0 0 -8 0 360/N 144/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 72/N 38/N 0 114/N 40/N 120/N 0 0 0 -8 0 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 6/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 76/N 0 40/N 0 -4 0 216/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 38/N 114/N 40/N 80/N 0 -4 540/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6/N*1i 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 40/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    4*1i 8/N*1i 0 8/N*1i 20/N*1i 12/N*1i 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -2 0 0 144/N 0 0 0 36/N 0 224/N 0 80/N 0 228/N 76/N 76/N 480/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 -2*1i 0 0 0 0 0 0 4/N*1i 12/N*1i 0 0 0 0 0 0 0 0;
    2*1i 4/N*1i 0 8/N*1i 0 6/N*1i 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -2 -2 -2 36/N 108/N 0 0 36/N 4/N -8/N 152/N 28/N 40/N 114/N 76/N 38/N -24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i -2*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0;
    2*1i 4/N*1i 0 8/N*1i 0 6/N*1i 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -2 -2 -2 36/N 0 108/N 0 36/N -4/N -8/N 152/N 36/N 40/N 114/N 76/N 30/N -24/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i -2*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0;
    2*1i 4/N*1i 0 8/N*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -2 -4 0 72/N 72/N 0 36/N 0 160/N -8/N 72/N -80/N -12/N 32/N 84/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 -4/N*1i -12/N*1i 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -4 -6 0 36/N 36/N 216/N 0 72/N -8/N 64/N 0 72/N -12/N -44/N -84/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0;
    4/N*1i 4*1i 0 0 0 0 0 16/N*1i 0 8/N*1i 16/N*1i 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 156/N 42/N 42/N -4/N 0 0 -2 -2 0 0 0 0 0 0 0 0 0 0 0 216/N 0 36/N 0 152/N 80/N 360/N 76/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 -4*1i 0 0 0 0 0 0 0 0 0 12/N*1i 0 0 0 0 0;
    2/N*1i 2*1i 0 0 0 0 0 8/N*1i 0 8/N*1i 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N 118/N -2/N 78/N 40/N -2 0 -4 -6 0 0 0 0 0 0 0 0 0 0 72/N 144/N 36/N 4/N 72/N 28/N -18/N 42/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2*1i 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0;
    2/N*1i 2*1i 0 0 0 0 0 8/N*1i 0 8/N*1i 0 6/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N -2/N 118/N 78/N 40/N -2 -4 0 -6 0 0 0 0 0 0 0 0 0 0 72/N 144/N 36/N -4/N 72/N 36/N -18/N 34/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i -2*1i 0 0 0 0 0 0 0 0 0 -6/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N -2/N 0 80/N 0 -2 -2 0 0 0 0 0 0 0 0 0 0 0 0 72/N 0 36/N -8/N 0 0 -44/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    2/N*1i 0 12*1i 0 0 0 0 0 20/N*1i 0 0 0 32/N*1i 36/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 38/N 40/N 40/N 38/N 0 0 0 0 0 -4 -4 0 0 0 0 0 0 0 0 0 0 144/N 0 0 0 0 0 224/N 0 228/N 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 2/N*1i 0 -4*1i 0 -4*1i 0 0 0 0 0 0 0 0 0 4/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 4/N -4/N 4/N 76/N 0 0 0 0 -4 -4 0 0 0 0 0 0 0 0 0 0 36/N 108/N 0 0 0 0 -16/N -4/N -132/N 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -8/N*1i 0 0;
    4/N*1i 0 0 4*1i 0 0 0 0 0 8/N*1i 0 0 0 0 16/N*1i 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 112/N -2/N -2/N 80/N 0 0 0 0 0 0 0 -2 -4 0 0 0 0 0 0 36/N 0 0 0 72/N 0 0 0 36/N 0 0 480/N 80/N 0 0 0 0 0 0 0 0 6/N*1i 0 2/N*1i 0 0 0 0 -10*1i 0 0 0 0 0 0 0 0 0 0 0;
    2/N*1i 0 0 2*1i 0 0 0 0 0 4/N*1i 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 74/N 74/N 0 40/N 0 0 0 0 0 0 -4 -6 0 0 0 0 0 0 0 36/N 0 0 72/N 0 0 0 36/N 0 0 -48/N -8/N 0 0 0 0 0 0 0 0 0 6/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    4/N*1i 0 0 0 10*1i 0 0 0 0 0 16/N*1i 0 0 0 0 24/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 76/N 30/N 38/N 76/N 4/N 0 0 0 0 0 0 0 0 -6 -4 0 0 0 0 0 0 36/N 0 0 108/N 0 0 0 4/N 0 0 152/N 0 0 0 0 0 0 0 0 0 4/N*1i 4/N*1i 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 38/N 38/N -76/N 76/N 0 0 0 0 0 0 0 0 -4 -2 0 0 0 0 0 0 0 36/N 0 36/N 0 0 0 0 0 0 -96/N 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0;
    4/N*1i 0 0 0 0 4*1i 0 0 0 0 0 8/N*1i 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 72/N 40/N 40/N 0 0 0 0 0 0 0 0 0 0 0 0 -6 0 0 0 0 0 0 0 72/N 0 108/N 0 0 0 36/N 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 4/N*1i 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0 0;
    6/N*1i 0 0 0 0 0 8*1i 0 0 0 0 12/N*1i 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 76/N 76/N 76/N 32/N -32/N 0 0 0 0 0 0 0 0 0 0 0 -6 0 0 0 0 0 0 0 72/N 0 72/N 72/N 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i -2/N*1i 0 0 0 0 0 0 0 -6*1i 0 0 0 0 0 0 0 0;
    8/N*1i 0 0 0 0 6*1i 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 44/N 36/N 80/N -80/N 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 0 0 0 0 0 108/N 0 -4/N 108/N 0 144/N 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0 0 0 0;
    2/N*1i 0 0 0 0 0 2*1i 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 0 0 0 108/N 0 0 0 0 144/N 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0 0 0 0;
    0 4/N*1i 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 24/N*1i 0 8/N*1i 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 118/N 40/N 40/N 0 0 0 40/N 0 0 0 0 0 0 0 2 -4 0 0 0 0 0 0 0 0 0 0 0 216/N 36/N 240/N 80/N 0 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0 -10*1i 0 0 0 0 0 0 0;
    0 2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 12/N*1i 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 78/N 78/N 126/N 0 0 0 40/N 0 0 0 0 0 0 -4 -2 0 0 0 0 0 0 0 0 0 0 0 216/N 36/N -24/N -8/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 2/N*1i 12/N*1i 0 0 0 0 0 10*1i 0 0 0 0 0 0 0 0 32/N*1i 0 0 24/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 38/N 40/N 40/N 0 156/N 42/N 0 0 40/N 0 0 0 0 0 0 0 -2 -4 0 0 0 0 0 0 0 0 0 0 216/N 0 0 152/N 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 6/N -2/N 114/N -4/N 122/N 0 0 0 40/N 0 0 0 0 0 0 -4 2 0 0 0 0 0 0 0 0 0 0 72/N 0 0 -96/N -4/N 0 0 0 0 0 2/N*1i 0 0 4/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0;
    0 4/N*1i 0 4/N*1i 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 16/N*1i 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 72/N 40/N 40/N 0 0 0 76/N 84/N 0 0 114/N 0 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 0 108/N 0 216/N 0 36/N 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0;
    0 4/N*1i 0 0 10/N*1i 0 0 0 0 0 8*1i 0 0 0 0 0 0 0 0 24/N*1i 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 76/N 32/N 40/N 12/N 0 0 0 0 118/N 46/N 0 76/N 0 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 0 72/N 0 144/N 0 4/N 0 0 0 0 0 4/N*1i 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0;
    0 2/N*1i 0 0 0 2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 0 0 0 0 0 0 0 0 0 42/N 0 0 152/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 144/N 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 0;
    0 6/N*1i 0 0 0 0 8/N*1i 0 0 0 0 6*1i 0 0 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 44/N 36/N -120/N 0 0 0 0 0 0 0 84/N 114/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 144/N 72/N -4/N 0 0 0 0 0 -2/N*1i 0 0 0 0 0 -2/N*1i 0 0 0 -6*1i 0 0 0 0 0;
    0 0 12/N*1i 2/N*1i 0 0 0 0 0 0 0 0 8*1i 0 0 0 0 0 0 0 12/N*1i 0 48/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 112/N 4/N 38/N 42/N 0 0 0 38/N 0 0 0 0 0 0 0 0 0 0 -4 0 0 0 0 0 36/N 0 0 72/N 0 0 0 0 0 0 0 8/N*1i 0 4/N*1i 10/N*1i 0 0 6/N*1i 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 0 4/N 0 0 0 -4/N 0 0 0 0 0 0 0 0 0 0 -12 0 0 0 0 0 0 0 0 108/N 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 -4*1i 0 0 0 0;
    0 0 12/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 12*1i 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N -40/N 0 0 0 0 42/N 0 38/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 72/N 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 2/N*1i 0 0 0 0 0 0 -4*1i 0 0;
    0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 0 0 0 0 0 0 72/N 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 36/N 0 0 0 0 0 0 0 0 0 0 0 0 10/N*1i 0 0 2/N*1i 0 0 0 0 0 0 0 0;
    0 0 0 4/N*1i 10/N*1i 0 0 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 0 80/N -40/N 0 0 72/N 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 0 36/N 36/N 0 0 0 0 0 0 0 0 0 0 0 6/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 32/N*1i 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 76/N 84/N 0 0 40/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 432/N 36/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 2/N*1i 10/N*1i 0 0 0 0 0 0 0 0 8*1i 0 0 0 0 0 0 36/N*1i 0 16/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 38/N 42/N 118/N 46/N 0 40/N 0 0 40/N 0 0 0 0 0 0 0 0 0 0 0 216/N 0 0 0 0 0 0 0 0 0 0 0 10/N*1i 2/N*1i 0 0 0 0 0 0;
    0 0 0 0 0 0 0 2/N*1i 0 2/N*1i 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40/N 0 0 0 42/N 0 114/N 0 0 0 0 80/N 0 0 0 2 0 0 0 216/N 0 0 0 0 0 0 0 0 0 0 0 0 10/N*1i 0 2/N*1i 0 0 0 0 0;
    0 0 0 0 0 0 0 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 24/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N 0 0 0 0 84/N 0 76/N 0 0 0 0 40/N 0 0 0 2 0 0 0 108/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 10/N*1i 2/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 32/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80/N -40/N 42/N 0 0 38/N 84/N 0 114/N 0 40/N 0 0 0 0 2 0 0 108/N 0 0 0 0 0 0 0 0 0 0 0 0 6/N*1i 6/N*1i 0 0 4/N*1i 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N 0 4/N 0 -4/N 0 114/N 0 0 0 0 0 0 0 0 -10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 -6*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 2*1i 0 0 0 20/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 42/N 0 40/N 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 8/N*1i 0 0 0 0 0 0 6*1i 0 0 0 32/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 42/N 84/N 0 40/N 40/N 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 4*1i 0 8/N*1i 4/N*1i 4/N*1i 0 0 0 8/N*1i 0 8/N*1i -4/N*1i 12/N*1i 4/N*1i 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 -4 0 72/N 36/N 0 0 400/N 80/N 80/N 120/N 0 0 0 0 0 0 0 -14*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 2*1i 0 2*1i 0 4/N*1i 4/N*1i 12/N*1i 0 0 0 8/N*1i 0 8/N*1i 0 4/N*1i 4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4 6 -4 72/N 36/N 4/N 0 -40/N -8/N -8/N -12/N 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 2*1i 2*1i 0 0 4/N*1i 4/N*1i 4/N*1i 0 0 0 8/N*1i 8/N*1i 4/N*1i 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4 2 72/N 0 -4/N 72/N 0 0 72/N 108/N 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 2/N*1i 2/N*1i 0 0 4*1i 2*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0 16/N*1i 16/N*1i 0 0 8/N*1i 8/N*1i 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 76/N 84/N 84/N 4 0 0 0 0 0 0 0 360/N 36/N 108/N 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 -2/N*1i 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 0 0 0 8/N*1i 4/N*1i 0 0 0 0 8/N*1i 0 12/N*1i 0 0 0 0 0 0 0 0 0 0 38/N 38/N 4/N 0 4 0 0 0 0 0 0 0 72/N 0 0 12/N 72/N 0 14/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N 0 0 4 0 0 0 0 0 0 0 0 -8/N 132/N 8/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 38/N 0 0 0 4 0 0 0 0 0 36/N 0 4/N 0 36/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 4/N*1i 0 0 0 0 0 0 0 0 40/N 0 0 0 0 0 0 6 0 0 0 36/N 0 0 0 0 0 0 14/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 4/N*1i -4/N*1i 0 0 0 0 0 0 0 0 4*1i -2*1i 0 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 80/N 0 0 0 0 0 0 0 6 0 0 0 36/N 0 4/N 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i 4/N*1i 4/N*1i 4/N*1i -4/N*1i 0 0 0 0 0 0 0 0 0 0 6*1i 0 4*1i 0 0 0 0 0 8/N*1i 0 0 8/N*1i 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 0 80/N 8/N 80/N 0 0 0 0 0 0 6 0 0 0 108/N -4/N 0 72/N 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 2/N*1i 2/N*1i -2/N*1i 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 8*1i 0 0 0 0 0 4/N*1i 12/N*1i 0 0 0 0 16/N*1i 0 0 0 0 0 0 0 0 0 40/N 4/N 40/N 0 0 0 0 0 0 0 6 0 0 72/N 0 -4/N 0 0 14/N*1i;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 12/N*1i 0 0 4/N*1i 0 0 0 0 0 0 0 42/N 0 0 0 40/N 0 0 0 10 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i -6/N*1i 4/N*1i -2/N*1i 0 0 0 0 0 0 0 0 0 0 4*1i -2*1i 0 0 0 0 0 0 0 0 0 0 16/N*1i 0 0 8/N*1i 0 0 0 0 0 0 42/N 84/N 0 84/N 0 40/N 0 0 0 10 0 0 0 0 12/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 2/N*1i -6/N*1i 0 0 0 0 0 0 0 0 2/N*1i 8/N*1i 0 0 0 0 0 0 6*1i 2*1i 0 0 0 0 0 0 0 16/N*1i 8/N*1i 0 0 0 0 0 0 0 44/N 0 0 0 0 0 42/N 84/N 0 0 12 0 0 0 -4/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4*1i 0 0 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 -8/N 8/N 0 4/N -4/N 0 0 0 0 -4 0 0 108/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i 0 0 0 0 0 0 4/N 40/N 0 0 0 0 -4/N 0 0 0 0 0 0 72/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N*1i -2/N*1i 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 0 0 0 6*1i 0 0 0 0 0 0 8/N*1i 0 0 0 0 0 0 0 44/N 4/N 44/N 0 0 42/N 0 0 0 0 0 0 12 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4/N -4/N 42/N 76/N 0 2 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 0 4/N*1i 0 0 0 0 4/N*1i 4/N*1i 0 0 0 0 0 0 0 0 14];
end
