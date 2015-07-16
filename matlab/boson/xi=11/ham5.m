function f = ham5(N)
f=[
    110-100/N 176/N 264/N 0 0 0 0 2*1i 0 0 0 0 0 0 0 0 0 0 0 0 0;
    110/N 88-100/N 0 264/N 176/N 0 0 0 2/N*1i 2*1i 0 0 0 0 0 0 0 0 0 0 0;
    220/N 0 110-100/N 44/N 176/N 0 0 -2/N*1i -2/N*1i 0 0 2*1i 0 0 0 0 0 0 0 0 0;
    0 88/N 0 66-100/N 0 264/N 0 0 0 0 0 0 0 2*1i 0 0 0 0 0 0 0;
    0 88/N 66/N 0 88-100/N 132/N 0 0 0 -2/N*1i 0 -2/N*1i 0 0 0 0 0 0 0 0 0;
    0 0 0 66/N 0 44-100/N 440/N 0 0 0 0 0 0 -2/N*1i 0 0 0 0 0 0 0;
    0 0 0 0 0 0 -100/N 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -10*1i -16/N*1i -24/N*1i 0 0 0 0 102-100/N 4 132/N 4/N 168/N 0 0 0 0 0 4*1i 0 0 0;
    0 0 0 0 0 0 0 4 98-100/N 44/N 0 -64/N 4/N 0 0 0 0 -2*1i 0 0 0;
    -10/N*1i -8*1i 0 -24/N*1i -16/N*1i 0 0 62/N 14/N 80-100/N 0 0 0 176/N 4/N 0 0 -2/N*1i 0 0 0;
    0 0 0 0 0 0 0 4/N 0 0 88-100/N 0 0 0 132/N 176/N 0 0 4*1i 0 0;
    -10/N*1i 0 -6*1i 0 -16/N*1i 0 0 80/N -40/N 0 0 98-100/N 0 44/N 0 4/N 0 -6/N*1i 0 0 0;
    0 0 0 0 0 0 0 0 4/N 0 0 0 110-100/N 0 44/N -88/N 0 0 0 6*1i 0;
    0 -8/N*1i 0 -6*1i 0 -24/N*1i 0 0 0 36/N 0 0 0 54-100/N 0 0 4/N 0 0 0 0;
    0 0 0 0 0 0 0 0 0 4/N 66/N 0 0 0 66-100/N 0 176/N 0 2/N*1i 0 6*1i;
    0 0 0 0 0 0 0 0 0 0 88/N 4/N -22/N 0 0 88-100/N 44/N 0 -4/N*1i 6/N*1i 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 4/N 44/N 0 44-100/N 0 0 0 0;
    0 0 0 0 0 0 0 -4*1i 2*1i -8/N*1i 0 -8/N*1i 0 0 0 0 0 90-100/N 4/N -12/N 0;
    0 0 0 0 0 0 0 0 0 0 -4*1i 0 0 0 -8/N*1i -8/N*1i 0 4/N 80-100/N 0 132/N;
    0 0 0 0 0 0 0 0 0 0 0 0 -2*1i 0 0 0 0 -4/N 0 98-100/N 44/N;
    0 0 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -2*1i 0 -8/N*1i 0 18/N 0 54-100/N];
end
