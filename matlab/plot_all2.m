function plot_all2(minBit, maxBit, points)
    if nargin < 3
        % by default, pick 1500 points.
        points = 1500;
    end
    addpath('..');

    xis = [1,2,5,10];
    for i = 1 : size(xis, 2)
        for bit = minBit : maxBit
            if bit == 3
                statenumber = 5;
            else
                statenumber = 8;
            end
            file = strcat('nxi=', num2str(xis(i)),'M=', num2str(bit), '.pdf');
            fprintf('Plotting %s ...\n', file);
            plot_states2(xis(i), bit, statenumber, {'points', points, 'file', file, 'mode', 'sr'});
        end
    end
end