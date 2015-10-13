function plot_groundex(xis)
    x = [3:2:11];
    y = zeros(1, size(x, 2));
    cols = ['r', 'g', 'b', 'm', 'c', 'k', 'y'];
    p = zeros(1, size(xis, 2));
    legendtitle = cell(1, size(xis, 2)); 
    figure; cla;
    hold on;
    for i = 1 : size(xis, 2)
        for k = 1 : size(x, 2)
            state = lowest_energies(xis(i), x(k), (x(k)-1)/2, 1);
            y(k) = state(1, 1);
        end
        p(i) = plot(x, y, strcat('--', cols(i), '.'));
        legendtitle{i} = strcat('xi=', num2str(xis(i)));
    end
    legend(p(:),legendtitle{:}, 'Location', 'southwest');
    xlabel('M');
    ylabel('E');
    hold off;
end