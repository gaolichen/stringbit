function plot_allhigh(minBit, maxBit, points)
    if nargin < 3
        % by default, pick 1500 points.
        points = 1500;
    end
    addpath('..');
    
    xis = [0,-1,1,2,5,10];
    for i = 1 : size(xis, 2)
        for bit = minBit : maxBit
            statenumber = 5;
            file = strcat('xi=', num2str(xis(i)),'M=', num2str(bit), '_lr', '.pdf');
            fprintf('Plotting %s ...\n', file);
            plot_states(xis(i), bit, statenumber, {'points', points, 'file', file, 'mode', 'lr'});
        end
    end
end