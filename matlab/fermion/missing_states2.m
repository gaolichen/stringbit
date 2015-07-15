function missing_states2(bits, statenumber)
    figs = {};
    for n = 2 : 2 * bits - 1
        res = lowest_energies2(bits, n / 2, statenumber);
        X = linspace(n/2 - 0.1, n/2 + 0.1);
        for i = 1 : size(res, 1)
            Y = res(i, 1);
            if res(i, 3) < 1e-8 || abs(res(i, 2)) > 1e-8 
                figs = cat(2, figs, X, Y, 'b');
            else
                figs = cat(2, figs, X, Y, 'r');
            end
        end
    end
    
    figure;
    plot(figs{:});
    title(strcat('Missing Energy Level for M=', num2str(bits)));
    ylabel('E');
    xlabel('N');
    xlim([0, bits + 1]);
end