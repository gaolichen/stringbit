function mat = norm_matrix(bits, fNumber, N)
%norm_matrix loads norm matrix for given bits and fermion number
%   Detailed explanation goes here
    [len, rawdata] = norm_raw_data(bits, fNumber);
    mat = zeros(len);
    %display(rawdata)
    %display(rawdata(1, 2))
    order = floor((bits + 1) / 2);
    %display(order)
    count = 1;
    for i = 1 : len
        for j = i : len
            lowOrd = rawdata(count, 1);
            %display(lowOrd)
            if lowOrd == 0
                mat(i, j) = 0;
            else
                if lowOrd == 2
                    lowOrd = 0;
                end
                
                val = 0;
                for k = 1 : order
                    val = val + rawdata(count, k + 1) / N^(lowOrd + 2 * k - 2);
                end
                %display(val)
                mat(i, j) = val;
            end
            
            if i ~= j
                mat(j, i) = mat(i, j);
            end
            count = count + 1;
        end
    end
end

