function zero_states(root)
%UNTITLED3 Summary of this function goes here
    cd(root);
    fid=fopen('zero-states.txt','w');

    addpath('../..');
    addpath('..');
    
    for bit = 3 : 11
        fprintf('processing %d bits...\n', bit);
        fprintf(fid, '%d\n', bit);
        for n = 1 : bit
            mat = metric(bit, n);
            cnt = 0;
            for i = 1 : size(mat, 1)
                if abs(mat(i, i)) > 1e-8
                    cnt = cnt + 1;
                end
            end
            
            fprintf(fid, '%d ', cnt);
        end
        fprintf(fid, '\n');
    end
    
    fclose(fid);


end

