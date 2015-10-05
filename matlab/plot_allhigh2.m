function plot_allhigh2(minBit, maxBit, points)
    if nargin < 3
        % by default, pick 1500 points.
        points = 1500;
    end
    addpath('..');
    
    xis = [1,2,5,10];
    for i = 1 : size(xis, 2)
        for bit = minBit : maxBit
            statenumber = 5;
            file = strcat('nxi=', num2str(xis(i)),'M=', num2str(bit), '_lr', '.pdf');
            fprintf('Plotting %s ...\n', file);
            plot_states2(bit, statenumber, points, 'lr', xis(i), file);
        end
    end
end