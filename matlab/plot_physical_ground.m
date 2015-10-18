function plot_physical_ground(xis, Ns, filename, legentLocations)
    if nargin < 3
        filename = '';
    end
    if nargin < 4
        legentLocations = {};
    end
    
    x = (3:2:11);
    y = zeros(1, size(x, 2));
    row = size(xis, 2) / 2;
    cols = ['r', 'g', 'b', 'm', 'c'];
    p = zeros(1, size(Ns, 2));
    legendtitle = cell(1, size(Ns, 2)); 
    fig = figure; cla;
    for i = 1 : size(xis, 2)
        subplot(row, 2, i);
        hold on;
        for j = 1 : size(Ns, 2)
            for k = 1 : size(x, 2)
                y(k) = real(physical_states(xis(i), x(k), Ns(j), 1));
            end
            p(j) = plot(x, y, strcat('--', cols(j), '.'));
            legendtitle{j} = strcat('N=', num2str(Ns(j)));
        end
        if size(legentLocations, 2) > 0
            legend(p(:),legendtitle{:}, 'Location', legentLocations{i});
        else
            legend(p(:),legendtitle{:}, 'Location', 'southwest');
        end
        xlabel('M');
        ylabel('E');
        
        if xis(i) == 0
            texTitle = '$$H = H_{0}$$';
        elseif xis(i) == 1
            texTitle = strcat('$$H = H_{0} + \Delta H$$');
        elseif xis(i) == -1
            texTitle = strcat('$$H = H_{0} - \Delta H$$');
        elseif xis(i) > 0
            texTitle = strcat('$$H = H_{0} +', num2str(xis(i)), '\Delta H$$');
        else
            texTitle = strcat('$$H = H_{0} ', num2str(xis(i)), '\Delta H$$');
        end
        
        title(texTitle, 'interpreter', 'latex');
        ax1 = gca;
        set(ax1, 'XTick', (3:2:11));
        hold off;
    end
    
    if ~strcmp(filename, '')
        saveas(fig, filename);
        close(fig);
    end
end