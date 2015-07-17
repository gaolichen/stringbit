function f = pham4(N)
f=[
    -48+40/N 0 0 0 0 -2*1i 0 0 0 0;
    48/N -36+40/N 0 0 0 0 -2*1i 0 0 0;
    0 0 -48+40/N 0 0 2/N*1i 0 0 0 0;
    0 36/N 48/N -24+40/N 0 0 2/N*1i 0 0 0;
    0 0 0 24/N 40/N 0 0 0 0 0;
    8*1i 12/N*1i 16/N*1i 0 0 -40+40/N 0 -4/N 0 0;
    8/N*1i 6*1i 0 16/N*1i 0 32/N -24+40/N 0 -4/N 0;
    0 0 0 0 0 -4/N 0 -36+40/N 0 -6*1i;
    0 0 0 0 0 0 -4/N 24/N -24+40/N 0;
    0 0 0 0 0 0 0 2*1i 4/N*1i -24+40/N];
end