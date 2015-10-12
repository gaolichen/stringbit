function plot_ground(xis, Ns)
    x = [3:2:11];
    y = zeros(1, size(x, 2));
    row = size(xis, 2) / 2;
    cols = ['r', 'g', 'b', 'm', 'c'];
    p = zeros(1, size(Ns, 2));
    legendtitle = cell(1, size(Ns, 2)); 
    figure; cla;
    for i = 1 : size(xis, 2)
        subplot(row, 2, i);
        hold on;
        for j = 1 : size(Ns, 2)
            for k = 1 : size(x, 2)
                state = lowest_energies(xis(i), x(k), Ns(j), 1);
                y(k) = state(1, 1);
            end
            p(j) = plot(x, y, strcat('--', cols(j), '.'));
            legendtitle{j} = strcat('N=', num2str(Ns(j)));
        end
        legend(p(:),legendtitle{:}, 'Location', 'southwest');
        xlabel('M');
        ylabel('E');
        title(strcat('$$\xi$$=', num2str(xis(i))), 'interpreter', 'latex');
        hold off;
    end
end