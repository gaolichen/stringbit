function norm_matrix_rank_all(maxBits)
%norm_matrix_rank_all output all rank data of norm matrix to a file.
%   Detailed explanation goes here

    fid=fopen('ranks_all.txt','w');
    
    for bit = 3 : maxBits
        fprintf('processing %d bits...\n', bit);
        for f = 0 : floor((bit-1)/2)
            fprintf(fid, '%d %d\n', bit, 2*f);
            for n = 1 : bit
                fprintf(fid, '%d ', rank(norm_matrix(bit, 2*f, n)));
            end
            fprintf(fid, '\n');
        end
    end
    
    fclose(fid);
end

