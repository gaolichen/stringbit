function plot_excitation(xis, Ns, filename, legentLocations)
%plot_excitation plot the excitation energy between ground physical states
% and first excitation states.
%
%   Detailed explanation goes here

    if nargin < 3
        filename = '';
    end
    if nargin < 4
        legentLocations = {};
    end
    
    x = (3:2:11);
    y = zeros(1, size(x, 2));
    if size(xis, 2) > 1
        row = size(xis, 2) / 2;
    else
        row = 0;
    end
    cols = ['r', 'g', 'b', 'm', 'k', 'c'];
    p = zeros(1, size(Ns, 2));
    legendtitle = cell(1, size(Ns, 2)); 
    fig = figure; cla;
    for i = 1 : size(xis, 2)
        if row > 0
            subplot(row, 2, i);
        end
        hold on;
        for j = 1 : size(Ns, 2)
            num = 0;
            for k = 1 : size(x, 2)
                states = real(physical_states(xis(i), x(k), Ns(j), 2));
                %fprintf('states size=%d\n', length(states));
                %disp(states);
                if length(states) < 2
                    break;
                else
                    y(k) = x(k)*(real(states(2) - states(1)));
                    num = num + 1;
                end
            end
            p(j) = plot(x(1:num), y(1:num), strcat('--', cols(j), '.'));
            legendtitle{j} = strcat('N=', num2str(Ns(j)));
        end
        if size(legentLocations, 2) >= i
            legend(p(:),legendtitle{:}, 'Location', legentLocations{i});
        elseif size(legentLocations, 2) == 0
            legend(p(:),legendtitle{:}, 'Location', 'southwest');
        end
        xlabel('M');
        ylabel('$$M\times (E_{1} - E{0})$$','interpreter', 'latex');
        
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
        
        if row == 0
            texTitle = strcat({'$$M\times (E_{1} - E_{0})$$ for '}, texTitle);
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

