function plot_all2(root, minBit, maxBit)
    cd(root);
    addpath('../..');
    addpath('..');
    
    folders ={'xi=1' 'xi=2' 'xi=5' 'xi=10'};
    points = 1000;
    for i = 1 : size(folders, 2)
        cd(folders{i});
        for bit = minBit : maxBit
            if bit == 3
                statenumber = 5;
            else
                statenumber = 8;
            end
            file = strcat('../', folders{i},'pM=', num2str(bit), '.pdf');
            fprintf('Plotting %s ...\n', file);
            plot_states2(bit, statenumber, points, file);
        end
        cd('..');
    end
end