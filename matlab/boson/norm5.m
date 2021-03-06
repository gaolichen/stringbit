function f = norm5(N)
A=[
    5+75/N^2+40/N^4, 20/N+100/N^3, 30/N+90/N^3, 60/N^2+60/N^4, 80/N^2+40/N^4, 120/N^3, 120/N^4;
    20/N+100/N^3, 4+68/N^2+48/N^4, 72/N^2+48/N^4, 24/N+96/N^3, 16/N+104/N^3, 72/N^2+48/N^4, 120/N^3;
    30/N+90/N^3, 72/N^2+48/N^4, 6+78/N^2+36/N^4, 6/N+114/N^3, 24/N+96/N^3, 48/N^2+72/N^4, 120/N^3;
    60/N^2+60/N^4, 24/N+96/N^3, 6/N+114/N^3, 6+78/N^2+36/N^4, 72/N^2+48/N^4, 36/N+84/N^3, 120/N^2;
    80/N^2+40/N^4, 16/N+104/N^3, 24/N+96/N^3, 72/N^2+48/N^4, 8+48/N^2+64/N^4, 24/N+96/N^3, 120/N^2;
    120/N^3, 72/N^2+48/N^4, 48/N^2+72/N^4, 36/N+84/N^3, 24/N+96/N^3, 12+108/N^2, 120/N;
    120/N^4, 120/N^3, 120/N^3, 120/N^2, 120/N^2, 120/N, 120
];

B=[
    1+2/N^2-3/N^4, -1/N^2+1/N^4, 3/N-3/N^3, 0, 4/N-4/N^3, 0, 6/N^2-6/N^4, 0, 0, 0;
    -1/N^2+1/N^4, 1-2/N^2+1/N^4, 1/N-1/N^3, 0, -2/N+2/N^3, 0, 2/N^2-2/N^4, 0, 0, 0;
    3/N-3/N^3, 1/N-1/N^3, 1+3/N^2-4/N^4, 0, 4/N^2-4/N^4, 0, 4/N-4/N^3, 0, 0, 0;
    0, 0, 0, 1+2/N^2-3/N^4, 0, -1/N^2+1/N^4, 0, 3/N-3/N^3, 4/N-4/N^3, 6/N^2-6/N^4;
    4/N-4/N^3, -2/N+2/N^3, 4/N^2-4/N^4, 0, 2+2/N^2-4/N^4, 0, 2/N-2/N^3, 0, 0, 0;
    0, 0, 0, -1/N^2+1/N^4, 0, 1-2/N^2+1/N^4, 0, 1/N-1/N^3, -2/N+2/N^3, 2/N^2-2/N^4;
    6/N^2-6/N^4, 2/N^2-2/N^4, 4/N-4/N^3, 0, 2/N-2/N^3, 0, 2-2/N^2, 0, 0, 0;
    0, 0, 0, 3/N-3/N^3, 0, 1/N-1/N^3, 0, 1+3/N^2-4/N^4, 4/N^2-4/N^4, 4/N-4/N^3;
    0, 0, 0, 4/N-4/N^3, 0, -2/N+2/N^3, 0, 4/N^2-4/N^4, 2+2/N^2-4/N^4, 2/N-2/N^3;
    0, 0, 0, 6/N^2-6/N^4, 0, 2/N^2-2/N^4, 0, 4/N-4/N^3, 2/N-2/N^3, 2-2/N^2
];

C=[
    1-5/N^2+4/N^4, 0, 0, 0;
    0, 1-2/N^2+1/N^4, 3/N^2-3/N^4, 3/N-3/N^3;
    0, 3/N^2-3/N^4, 3-12/N^2+9/N^4, 3/N-3/N^3;
    0, 3/N-3/N^3, 3/N-3/N^3, 3-3/N^2
];

f=blkdiag(A, B, C);
end
