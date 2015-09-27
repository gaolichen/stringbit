function plot_allhigh2(root, minBit, maxBit, points)
    if nargin < 4
        % by default, pick 1500 points.
        points = 1500;
    end
    cd(root);
    addpath('../..');
    addpath('..');
    
    folders ={'xi=1' 'xi=2' 'xi=5' 'xi=10'};
    xis = [1,2,5,10];
    for i = 1 : size(folders, 2)
        cd(folders{i});
        for bit = minBit : maxBit
            statenumber = 5;
            file = strcat('../', folders{i},'pM=', num2str(bit), '_lr', '.pdf');
            fprintf('Plotting %s ...\n', file);
            plot_states2(bit, statenumber, points, 'lr', xis(i), file);
        end
        cd('..');
    end
end