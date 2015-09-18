function plot_all(root, minBit, maxBit, points)
    if nargin < 4
        % by default, pick 1500 points.
        points = 1500;
    end

    cd(root);
    addpath('../..');
    addpath('..');
    
    folders ={'xi=0' 'xi=n1' 'xi=1' 'xi=2' 'xi=5' 'xi=10'};
    xis = [0,-1,1,2,5,10];
    for i = 1 : size(folders, 2)
        cd(folders{i});
        for bit = minBit : maxBit
            if bit == 3
                statenumber = 5;
            else
                statenumber = 8;
            end
            file = strcat('../', folders{i},'M=', num2str(bit), '.pdf');
            fprintf('Plotting %s ...\n', file);
            plot_states(bit, statenumber, points, xis(i), file);
        end
        cd('..');
    end
end