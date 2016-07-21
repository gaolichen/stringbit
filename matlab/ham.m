function f = ham(bits, n)
    file1 = strcat(num2str(bits),'h0Ham-re-0.dat');
    file2 = strcat(num2str(bits),'h0Ham-re-1.dat');
    file3 = strcat(num2str(bits),'h0Ham-im-0.dat');
    file4 = strcat(num2str(bits),'h0Ham-im-1.dat');
    
    mat1 = load(file1);
    mat2 = load(file2);
    mat3 = load(file3);
    mat4 = load(file4);
    f = spconvert(mat1) + spconvert(mat2)/n + i*spconvert(mat3) + i*spconvert(mat4)/n;
end