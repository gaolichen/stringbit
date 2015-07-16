function f = ham4(N)
f=[
    -32+40/N 0 0 0 0 2*1i 0 0 0 0;
    32/N -24+40/N 0 0 0 0 2*1i 0 0 0;
    0 0 -32+40/N 0 0 -2/N*1i 0 0 0 0;
    0 24/N 32/N -16+40/N 0 0 -2/N*1i 0 0 0;
    0 0 0 16/N 40/N 0 0 0 0 0;
    -8*1i -12/N*1i -16/N*1i 0 0 -40+40/N 0 4/N 0 0;
    -8/N*1i -6*1i 0 -16/N*1i 0 8/N -36+40/N 0 4/N 0;
    0 0 0 0 0 4/N 0 -24+40/N 0 6*1i;
    0 0 0 0 0 0 4/N 16/N -16+40/N 0;
    0 0 0 0 0 0 0 -2*1i -4/N*1i -36+40/N];
end
