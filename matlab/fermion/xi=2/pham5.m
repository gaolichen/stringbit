function f = pham5(N)
f=[
    -20+10/N 0 0 0 0 0 0 0 0 0 0 0 -4*1i -2*1i 0 0 0 0 0 0 0;
    16/N -16+10/N 0 0 0 0 0 0 0 0 0 0 0 -4/N*1i -4*1i 0 0 0 0 0 0;
    4/N 0 -16+10/N 0 0 0 0 0 0 0 0 0 2/N*1i 0 0 -2*1i 0 0 0 0 0;
    0 0 0 -20+10/N 0 0 0 0 0 0 0 0 2/N*1i 2/N*1i 0 0 -6*1i 0 0 0 0;
    0 0 0 0 -20+10/N 0 0 0 0 0 0 0 0 4/N*1i 0 0 0 -2*1i 0 0 0;
    0 12/N 0 8/N 0 -12+10/N 0 0 0 0 0 0 0 0 -2/N*1i 0 0 0 -6*1i 0 0;
    0 4/N 16/N 0 8/N 0 -12+10/N 0 0 0 0 0 0 0 2/N*1i 0 0 0 0 -2*1i 0;
    0 0 0 8/N 12/N 0 0 -16+10/N 0 0 0 0 0 0 4/N*1i 0 0 2/N*1i 0 0 0;
    0 0 0 4/N 0 0 0 0 -16+10/N 0 0 0 0 0 0 2/N*1i 6/N*1i 0 0 0 0;
    0 0 0 0 0 8/N 0 8/N 0 -8+10/N 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 4/N 12/N 8/N 16/N 0 -8+10/N 0 0 0 0 0 0 0 6/N*1i 2/N*1i 0;
    0 0 0 0 0 0 0 0 0 8/N 8/N 10/N 0 0 0 0 0 0 0 0 0;
    4*1i 8/N*1i 0 8/N*1i 12/N*1i 0 0 0 0 0 0 0 -8+10/N -4 0 0 24/N 8/N 0 0 -10*1i;
    2*1i 4/N*1i 0 8/N*1i 0 0 0 0 0 0 0 0 -4 -12+10/N 0 0 -24/N -8/N 0 0 0;
    4/N*1i 4*1i 0 0 0 16/N*1i 0 8/N*1i 0 0 0 0 8/N 16/N -8+10/N 0 0 0 0 0 0;
    2/N*1i 0 8*1i 0 0 0 12/N*1i 0 16/N*1i 0 0 0 4/N 8/N 0 -8+10/N 0 0 0 0 10/N*1i;
    2/N*1i 0 0 2*1i 0 0 0 4/N*1i 0 0 0 0 4/N 0 0 0 -8+10/N 0 0 0 10/N*1i;
    4/N*1i 0 0 0 6*1i 0 0 8/N*1i 0 0 0 0 8/N 0 0 0 0 -8+10/N 0 0 0;
    0 2/N*1i 0 0 0 2*1i 0 0 0 12/N*1i 0 0 0 0 8/N 0 8/N 0 10/N 0 0;
    0 2/N*1i 8/N*1i 0 0 0 6*1i 0 0 0 16/N*1i 0 0 0 8/N 16/N 0 8/N 0 10/N 0;
    0 0 0 0 0 0 0 0 0 0 0 0 2*1i 0 4/N*1i 0 0 4/N*1i 0 0 10/N];
end
