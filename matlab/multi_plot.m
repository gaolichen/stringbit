function multi_plot(xis, bits, statenumber, titles, args)
    if nargin < 5
        args = {};
    end
    
    filename = '';
    if size(args, 2) > 0
        filename = GetArg(args, 'file', filename);
    end
    
    args = [args {'subplot', 1}];
    n = size(xis, 2) * size(bits, 2);
    row = n / 2;
    cnt = 0;
    fig = figure; cla;
    for i = 1 : size(xis, 2)
        for j = 1 : size(bits, 2)
            cnt = cnt + 1;
            if mod(cnt, 2) == 1
                h1 = subplot(row, 2, cnt);
            else
                h2 = subplot(row, 2, cnt);
            end
            
            plot_states(xis(i), bits(j), statenumber, args);
            
            title(titles{cnt});
            if mod(cnt, 2) == 1
                ylabel('E');
            end
            
            if n - cnt <= 1
                xlabel('1/N');
            end
            
            if mod(cnt, 2) == 0
                p1 = get(h1, 'pos');
                p2 = get(h2, 'pos');
                dx = (p2(1) - p1(1) - p1(3)) / 3;
                p1(3) = p1(3) + dx;
                p2(1) = p2(1) - dx;
                p2(3) = p2(3) + dx;
                set(h1, 'pos', p1);
                set(h2, 'pos', p2);
            end
        end
    end
    
    if ~strcmp(filename, '')
        saveas(fig, filename);
        close(fig);
    end
end

function f = GetArg(argList, name, defaultValue)
    for i = 1 : 2 : size(argList, 2)
        if strcmp(argList{i}, name)
            f = argList{i+1};
            return
        end
    end
    f = defaultValue;
end