function res=find_match(mat, n)
    mat = sortrows(mat, 1);
    res = zeros(1, n);
    visited = zeros(1, n);
    for i = 1 : size(mat, 1)
        a = int32(mat(i, 2));
        b = int32(mat(i, 3));
        if res(a) == 0 && visited(b) == 0
            res(a) = b;
            visited(b) = 1;
        end
    end    
end