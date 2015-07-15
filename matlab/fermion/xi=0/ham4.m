function f = ham4(N)
f=[
    8/N 0 0 0 0 0 0 4*1i 0 0;
    0 8/N 0 0 0 0 0 2/N*1i 6*1i 0;
    0 0 8/N 0 0 0 0 -2/N*1i 0 2*1i;
    0 0 0 8/N 0 0 0 -4/N*1i 0 0;
    0 0 0 0 8/N 0 0 0 0 0;
    0 0 0 0 0 8/N 0 0 -6/N*1i -2/N*1i;
    0 0 0 0 0 0 8/N 0 0 0;
    -4*1i -8/N*1i 0 -8/N*1i 0 0 0 -8+8/N 0 0;
    -2/N*1i -2*1i 0 0 -8/N*1i 0 0 -4/N -12+8/N 0;
    -2/N*1i 0 -6*1i 0 0 -8/N*1i 0 -4/N 0 -12+8/N];
end
