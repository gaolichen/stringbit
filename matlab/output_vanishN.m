function output_vanishN(root)
    cd(root);
    fid=fopen('ground_vanish_N.txt','w');

    addpath('../..');
    addpath('..');

    folders ={'xi=0' 'xi=n1' 'xi=1' 'xi=2' 'xi=5' 'xi=10'};
    %folders ={'xi=0'};

    for i = 1 : size(folders, 2)
        cd(folders{i});
        
        fprintf('processing %s ...\n', folders{i});
        xi = folders{i}(4:length(folders{i}));
        xi = strrep(xi, 'n', '-');

        fprintf(fid, '%s\n', xi);
        fprintf('ground_vanish_N\n');
        fprintf(fid, 'ground_vanish_N\n');
        res = ground_vanish_N(11);
        for j = 1 : size(res, 1)
            fprintf(fid, '%2.0f %4.1f %f %f %f\n', res(j, 1), res(j, 2), res(j, 3), res(j, 4), res(j, 5));
        end

        if i > 2
            fprintf('ground_vanish_N2\n');
            fprintf(fid, 'ground_vanish_N2\n');
            res = ground_vanish_N2(11);
            for j = 1 : size(res, 1)
                fprintf(fid, '%2.0f %4.1f %f %f %f\n', res(j, 1), res(j, 2), res(j, 3), res(j, 4), res(j, 5));
            end
        end

        cd('..');
    end

    fclose(fid);
end
