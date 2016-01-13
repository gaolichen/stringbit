function multi_plot2(xis, bits, statenumber, subtitles, args, maxx)
    if nargin < 5
        args = {};
    end
    if nargin < 6
        maxx = {};
    end
    
    filename = '';
    texTitle = '';
    pos = [];
    col = 2;
    if size(args, 2) > 0
        filename = GetArg(args, 'file', filename);
        texTitle = GetArg(args, 'title', texTitle);
        pos = GetArg(args, 'position', pos);
        col = GetArg(args, 'column', col);
    end
    
    args = [args {'subplot', 1}];
    n = size(xis, 2) * size(bits, 2);
    row = n / col;
    cnt = 0;
    fig = figure; cla;
    if size(pos, 2) == 4
        set(fig, 'Position',pos);
    end
    
    % If more than two row, increase the size of the figure.
    if row > 2
        ppos = get(fig, 'PaperPosition');
        ppos(4) = (ppos(4) - ppos(2)) * row / 2 + ppos(2);
        set(fig, 'PaperPosition', ppos);
    end

    for i = 1 : size(xis, 2)
        for j = 1 : size(bits, 2)
            cnt = cnt + 1;
            if mod(cnt, 2) == 1
                h1 = subplot(row, 2, cnt);
            else
                h2 = subplot(row, 2, cnt);
            end
            
            if size(maxx, 2) > 0
                plot_states2(xis(i), bits(j), statenumber, [args {'maxx', maxx(cnt)}]);
            else
                plot_states2(xis(i), bits(j), statenumber, args);
            end
            
            title(subtitles{cnt}, 'interpreter', 'latex');
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
    
    if ~strcmp(texTitle, '')
        text(-1, 1, texTitle, 'interpreter', 'latex');
    end

    if ~strcmp(filename, '')
        saveas(fig, filename);
        close(fig);
    end
end

function f = GetArg(argList, name, defaultValue)
    for i = 1 : 2 : size(argList, 2)
        if strcmpi(argList{i}, name)
            f = argList{i+1};
            return
        end
    end
    f = defaultValue;
end