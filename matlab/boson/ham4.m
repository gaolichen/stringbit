function f = ham4(N)
f=[
    8 12/N 16/N 0 0 2*1i 0 0 0 0;
    8/N 6 0 16/N 0 0 2*1i 0 0 0;
    8/N 0 8 4/N 0 -2/N*1i 0 0 0 0;
    0 6/N 0 4 24/N 0 -2/N*1i 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    -8*1i -12/N*1i -16/N*1i 0 0 0 8/N 4/N 0 0;
    -8/N*1i -6*1i 0 -16/N*1i 0 -4/N -6 0 4/N 0;
    0 0 0 0 0 4/N 0 6 8/N 6*1i;
    0 0 0 0 0 0 4/N 4/N 4 0;
    0 0 0 0 0 0 0 -2*1i -4/N*1i -6];
end
