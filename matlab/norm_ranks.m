function norm_ranks(root)
%norm_ranks calculate ranks of all norm matrix.

    cd(root);
    fid=fopen('ranks.txt','w');

    addpath('../..');
    addpath('..');
    
    for bit = 3 : 11
        fprintf('processing %d bits...\n', bit);
        fprintf(fid, '%d\n', bit);
        for n = 1 : bit
            fprintf(fid, '%d ', rank(metric(bit, n)));
        end
        fprintf(fid, '\n');
    end
    
    fclose(fid);
end

