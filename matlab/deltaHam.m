function f = deltaHam(bits, n)
    file = strcat(num2str(bits),'deltaHam-re-1.dat');
    mat = load(file);
    f = spconvert(mat)/n;
end