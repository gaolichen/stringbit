function f = plot_states(xi, bits, statenumber, args)
    if nargin < 4
        args = {};
    end
    
    % by default, pick 100 points.
    points = 100;
    mode = 'sr';
    filename = '';
    issubplot = 0;
    maxX = 1.5;
        
    if size(args, 2) > 0
        points = GetArg(args, 'points', points);
        mode = GetArg(args, 'mode', mode);
        filename = GetArg(args, 'file', filename);
        issubplot = GetArg(args, 'subplot', issubplot);
        maxX = GetArg(args, 'maxx', maxX);
    end
        
    minN = 1/maxX;
    % X: the x-components
    X = 0:1/minN/points:1/minN;
    
    tot = length(X);
    % stored all the states.
    allstates = zeros(tot, statenumber, 3);
    % whether the current state is degenerate.
    isDegenerate = zeros(1, statenumber, 'int8');
    % pos1, pos2 : the index to calculate next point.
    pos1 = zeros(1, statenumber, 'int16');
    pos2 = zeros(1, statenumber, 'int16');
    % the distance matrix;
    dismat = zeros(statenumber * statenumber, 3);
    
    % midX: pick a point other than X(1) to start to avoid degenerate case.
    %midX = min(10, round((1/(2*bits)/maxX) * points));
    midX = round((1/(2*bits)/maxX) * points);
    newPlotStart = 0;
    
    % first go through midX to tot and connect lines.
    for curr = midX : tot
        n = 1/X(curr);
        states = lowest_energies(xi, bits, n, statenumber, mode);
        initDegenerate;
        
        if curr == midX
            allstates(midX,1:statenumber,:) = states(1:statenumber,:);
            initHighestE = states(statenumber, 1);
            initLowestE = states(1, 1);
            highestE = initHighestE;
            lowestE = initLowestE;
            initGap = initHighestE - initLowestE;
        else
            build_dismat;
            res = find_match(dismat, statenumber);
            for j = 1 : statenumber
                allstates(curr, j, :) = states(res(j), :);
            end
            
            highestE = max(states(statenumber, 1), highestE);
            lowestE = min(states(1, 1), lowestE);
            finalHighE = states(statenumber, 1);
            finalLowE = states(1, 1);
            currGap = highestE - lowestE;
            
            if currGap > 3.5 * initGap && newPlotStart == 0
                newPlotStart = curr;
                midHighE = states(statenumber, 1);
                midLowE = states(1, 1);
            end
        end
        
        nextPos;
    end
    
    fprintf('initGap=%d currGap=%d\n', initGap, currGap);
    if currGap <= 3.5 *initGap
        midHighE = finalHighE;
        midLowE = finalLowE;
    end
    
    % set pos1=midX+1 and pos2 = midX
    for i = 1 : statenumber
        pos1(i) = midX + 1;
        pos2(i) = midX;
    end
    
    % go through from midX -1 to 1.
    for curr = midX - 1 : -1 : 1
        if X(curr) == 0
            n = Inf;
        else
            n = 1/X(curr);
        end
        
        states = lowest_energies(xi, bits, n, statenumber, mode);
        initDegenerate;
        build_dismat;
        res = find_match(dismat, statenumber);
        for j = 1 : statenumber
            allstates(curr, j, :) = states(res(j), :);
        end
        
        nextPos;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % nested function.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set isDegenerate array.
    function initDegenerate
        ii = 1;
        while ii < statenumber
            if IsEqual(states(ii, 1), states(ii + 1, 1)) && IsEqual(states(ii, 2), states(ii + 1, 2))
                isDegenerate(ii) = 1;
                isDegenerate(ii + 1) = 1;
                ii = ii + 2;
            else
                isDegenerate(ii) = 0;
                ii = ii + 1;
            end
        end
    end

    % return the next possible point for given state
    function f = nextPoint(p1, p2, stateid)
        if p2 == p1
            f = allstates(p1, stateid, :);
        else
            f = allstates(p2, stateid, :) + (allstates(p2, stateid, :) - allstates(p1, stateid, :)) * double(curr - p2) / double(p2 - p1);
        end
    end

    % return next pos1, pos2 for all states.
    function nextPos
        for ii = 1 : statenumber
            if pos1(ii) == 0
                pos1(ii) = curr;
                pos2(ii) = curr;
            elseif ~isDegenerate(ii)
                pos1(ii) = pos2(ii);
                pos2(ii) = curr;
            end
        end
    end

    % build the distance matrix for states.
    function build_dismat
        for ii = 1 : statenumber
            next = nextPoint(pos1(ii), pos2(ii), ii);
            %fprintf('ii = %d pos1=%d pos2=%d next size=%d\n', ii, pos1(ii), pos2(ii), size(next, 1));
            for jj = 1 : statenumber
                %fprintf('jj size=%d\n',  size(squeeze(states(jj,:)),2));
                if isDegenerate(jj)
                    %fprintf('curr = %d ii=%d jj=%d kk=%d \n', curr, ii, jj, (ii-1)*statenumber + jj);
                    dismat((ii-1)*statenumber + jj, 1) = norm([next(1)-states(jj,1), next(2)-states(jj,2)]);
                else
                    dismat((ii-1)*statenumber + jj, 1) = norm([next(1)-states(jj,1), next(2)-states(jj,2), next(3)-states(jj,3)]);
                end
                
                dismat((ii-1)*statenumber + jj, 2) = ii;
                dismat((ii-1)*statenumber + jj, 3) = jj;
            end
        end
    end

    % break the lines into solid/dotted/dashed parts.
    lines(1, bits * statenumber) = PlotInfo();
    lineNumber = 0;
    for staid = 1  : statenumber
        % lineType:
        % 0 for unset; 
        % 1 for solid line; 
        % 2 for dash-dot (complex) line
        % 3 for dotted line
        lineType = 1;
        sx = 1;
        i = 1;
        while i <= tot
            if lineType == getLineType(staid, i)
                i = i + 1;
            else
                % save the current line
                makenewline;
                if lineType == 2
                    % for complex (dashed) line, we need to check
                    % duplicates.
                    for j = lineNumber - 1 : -1 : 1
                        if lines(j).LineType == 2 && lines(j).StateId ~= staid
                            % check if two lines cover each other.
                            if covers(j, lineNumber)
                                lines(j) = lines(j).MixState(staid);
                                lineNumber = lineNumber - 1;
                                break;
                            elseif covers(lineNumber, j)
                                sid = lines(j).StateId;
                                lines(j) = lines(lineNumber).MixState(sid);
                                lineNumber = lineNumber - 1;
                                break;
                            end
                        end
                    end
                end
                lineType = getLineType(staid, i);
                i = i + 1;
            end
        end
        
        % save last line.
        if sx < tot
            i = i - 1;
            makenewline;
        end
    end
    
    fprintf('total line segments: %d\n', lineNumber);

    % return if two line segments cover each other.
    function f = covers(l1, l2)
        a = lines(l1);
        b = lines(l2);
        if a.StartX <= b.StartX && a.EndX >= b.EndX
            f = isConjugate(a.StateId, b.StateId, b.StartX + 1) && isConjugate(a.StateId, b.StateId, b.EndX - 1);
        else
            f = 0;
        end
    end
    
    % make new line segments.
    function makenewline
        lineNumber = lineNumber + 1;  
        lines(lineNumber) = PlotInfo(sx, i, staid, lineType);
        sx = i;
    end
    
    % return if two points on two lines are complex conjugate to each other.
    function f = isConjugate(sid1, sid2, xpos)
        f = IsEqual(allstates(xpos, sid1, 1), allstates(xpos, sid2, 1)) && IsEqual(allstates(xpos, sid1, 2), -allstates(xpos, sid2, 2));
    end
    
    function f = getLineType(stateid, xpos)
        im = allstates(xpos, stateid, 2);
        nm = allstates(xpos, stateid, 3);
        
        if ~IsEqual(im, 0)
            f = 2;
        elseif IsNegative(nm)
            f = 3;
        else
            f = 1;
        end
    end

    if strcmp(mode, 'sr')
        plainTitle = strcat({'Lowest '},  num2str(statenumber), {' Energy Levels for '});
    else
        plainTitle = strcat({'Highest '},  num2str(statenumber), {' Energy Levels for '});
    end

    if xi > -10
        if xi == 0
            texTitle = '$$H = H_{0},';
        elseif xi == 1
            texTitle = strcat('$$H = H_{0} + \Delta H,');
        elseif xi == -1
            texTitle = strcat('$$H = H_{0} - \Delta H,');
        elseif xi > 0
            texTitle = strcat('$$H = H_{0} +', num2str(xi), '\Delta H,');
        else
            texTitle = strcat('$$H = H_{0} ', num2str(xi), '\Delta H,');
        end
        
        texTitle = strcat(texTitle, '\,M = ', num2str(bits), '$$');
    else
        texTitle = strcat('$$M = ', num2str(bits), '$$');
    end
    
    texTitle = strcat(plainTitle, texTitle);

    % plot results.
    if issubplot == 0
        fig = figure; cla;
    end
    
    if issubplot || currGap < initGap * 8
        doplot(1, tot, true);
    else
        fprintf('newPlotStart = %d\n', newPlotStart);
        h1 = subplot(1, 2, 1);
        doplot(1, newPlotStart, false);
        h2 = subplot(1, 2, 2);
        doplot(newPlotStart, tot, false);
        
        p1 = get(h1, 'pos');
        p2 = get(h2, 'pos');
        dx = (p2(1) - p1(1) - p1(3)) / 3;
        p1(3) = p1(3) + dx;
        p2(1) = p2(1) - dx;
        p2(3) = p2(3) + dx;
        set(h1, 'pos', p1);
        set(h2, 'pos', p2);
        axes('Position',[0 0 1 1],'Visible','off');
        text(0.35,0.95, texTitle, 'interpreter', 'latex');
    end
    
    if issubplot ==0 && ~strcmp(filename, '')
        saveas(fig, filename);
        close(fig);
    end
    
    function doplot(startX, endX, singleplot)
        p = zeros(1, 3);
        legendtitle = {'positive norm', 'zero norm', 'negative norm'};
        
        hold on;
        for ii = 1 : lineNumber
            info = lines(ii);
            from = max(startX, info.StartX);
            to = min(endX, info.EndX);
            if from < to
                p(info.LineType) = plot(X(from:to), allstates(from:to, info.StateId, 1), 'Color', info.Color, 'LineStyle', info.LineStyle, 'LineWidth', info.LineWidth);
            end
        end
        
        for ii = 3 : -1 : 1
            if p(ii) == 0
                p(ii) = [];
                legendtitle(ii) = [];
            end
        end
        
        % do not generate legend is it is called by multi_plot function.
        if issubplot == 0
            % place the legendLocation.
            if startX == 1
                highE1 = initHighestE;
                highE2 = midHighE;
                lowE1 = initLowestE;
                lowE2 = midLowE;
            else
                highE1 = midHighE;
                highE2 = finalHighE;
                lowE1 = midLowE;
                lowE2 = finalLowE;
            end

            if abs(highE1-highE2) > abs(lowE1-lowE2)
                if highE1 > highE2
                    legendLocation = 'northeast';
                else
                    legendLocation = 'northwest';
                end
            else
                if lowE1 > lowE2
                    legendLocation = 'southwest';
                else
                    legendLocation = 'southeast';
                end
            end
            
            legend(p(:),legendtitle{:}, 'Location', legendLocation);
        end
        
        ax1 = gca;
        
        % dont generate lable for multi_plot
        if issubplot == 0
            xlabel('1/N');
            if startX == 1
                ylabel('E');    
            else
                %set(ax1, 'YAxisLocation', 'right');
            end
        end
        
        if singleplot && ~issubplot
            title(texTitle, 'interpreter', 'latex');
            set(ax1, 'XTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]);
        else
        end
        
        xlim([X(startX), X(endX)]);
        ylim auto;
        hold off
    end
end

function f = IsEqual(a, b)
    f = abs(a - b) < 1e-5;
end
   
function f = IsNegative(a)
    f = (a < -1e-6);
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