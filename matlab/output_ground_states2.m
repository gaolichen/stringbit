function output_ground_states2(root)
    cd(root);
    fid=fopen('ground_states2.txt','w');

    addpath('../..');
    addpath('..');

    folders ={'xi=1' 'xi=2' 'xi=5' 'xi=10'};
    %folders ={'xi=0'};

    for i = 1 : size(folders, 2)
        cd(folders{i});
        xi = folders{i}(4:length(folders{i}));
        xi = strrep(xi, 'n', '-');
        fprintf('processing %s ...\n', folders{i});
        fprintf(fid, '%s\n', xi);

        for bit = 3 : 11
            [e, f0, f1] = get_ground_state2(bit, 10000);
            fprintf(fid, '%d %d %d\n', bit, size(f0, 1), size(f1, 1));
            fprintf(fid, '%f\t%f\n', real(e(1)), real(e(2)));

            fprintf(fid, '\n');
            for j = 1 : size(f0, 1)
                fprintf(fid, '%d\t%f\n', round(f0(j, 1)), f0(j, 3));
            end

            fprintf(fid, '\n');
            for j = 1 : size(f1, 1)
                fprintf(fid, '%d\t%f\n', round(f1(j, 1)), f1(j, 3));
            end

            fprintf(fid, '\n');
        end

        cd('..');
    end

    fclose(fid);
end