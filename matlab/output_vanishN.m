function output_vanishN(root)
    cd(root);
    fid=fopen('ground_vanish_N.txt','w');

    addpath('../..');
    addpath('..');

    folders ={'xi=0' 'xi=1' 'xi=5' 'xi=10' 'xi=n1' 'xi=n5' 'xi=n10'};
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
            fprintf(fid, '%2.0f %4.1f\n', res(j, 1), res(j, 2));
        end

        fprintf('ground_vanish_N2\n');
        fprintf(fid, 'ground_vanish_N2\n');
        res = ground_vanish_N2(11);
        for j = 1 : size(res, 1)
            fprintf(fid, '%2.0f %4.1f\n', res(j, 1), res(j, 2));
        end

        cd('..');
    end

    fclose(fid);
end
