function plot_states4(bits, statenumber, points, maxX, filename)
    if nargin < 3
        % by default, pick 100 points.
        points = 100;
    end
    if nargin < 4
        maxX = 1.5;
    end
    if nargin < 5
        filename = '';
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
    midX = round((1/(1.5*bits)/maxX) * points);
    newPlotStart = 0;
    
    % first go through midX to tot and connect lines.
    for curr = midX : tot
        n = 1/X(curr);
        states = lowest_energies(bits, n, statenumber);
        initDegenerate;
        
        if curr == midX
            allstates(midX,1:statenumber,:) = states(1:statenumber,:);
            highestE = states(statenumber, 1);
            lowestE = states(1, 1);
            initGap = highestE - lowestE;
        else
            build_dismat;
            res = find_match(dismat, statenumber);
            for j = 1 : statenumber
                allstates(curr, j, :) = states(res(j), :);
            end
            
            highestE = max(states(statenumber, 1), highestE);
            lowestE = min(states(1, 1), lowestE);
            currGap = highestE - lowestE;
            
            if currGap > 4 * initGap && newPlotStart == 0
                newPlotStart = curr;
            end
        end
        
        nextPos;
    end
    fprintf('initGap=%d currGap=%d\n', initGap, currGap);
    
    % set pos1 and pos2 = midX
    for i = 1 : statenumber
        pos1(i) = midX;
        pos2(i) = midX;
    end
    
    % go through from midX -1 to 1.
    for curr = midX - 1 : -1 : 1
        if X(curr) == 0
            n = Inf;
        else
            n = 1/X(curr);
        end
        
        states = lowest_energies(bits, n, statenumber);
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
    lines(1, 10 * statenumber) = PlotInfo();
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
                    for j = lineNumber : -1 : 1
                        if lines(j).LineType == 2 && abs(lines(j).StartX - i) <= 5 
                            % check if the eigenvalues are conjugate to
                            % each other.
                            sid = lines(j).StateId;
                            if isConjugate(sid, staid, i + 10)
                                % this is not a new complex line.
                                lines(j) = lines(j).MixState(staid);
                                i = lines(j).EndX;
                                sx = i;
                                lineType = getLineType(staid, i);
                                break;
                            end
                        end
                    end
                    
                    %fprintf('j=%d\n', j);
                    if j == 1
                        fprintf('staid=%d, i=%d, re=%f, im=%f\n', staid, i, allstates(i, staid, 1), allstates(i, staid, 2));
                    end
                end
                
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
%     for i = 1 : lineNumber
%         disp(lines(i));
%     end
    
    function makenewline
        lineNumber = lineNumber + 1;  
        lines(lineNumber) = PlotInfo(sx, i, staid, lineType);
        lineType = getLineType(staid, i);
        sx = i;
    end
    
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

    % plot results.
    fig = figure; cla;
    if currGap < initGap * 9
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
        text(0.45,0.95,strcat({'Lowest '},  num2str(statenumber), ' Energy States for M=', num2str(bits)));
    end
    
    if ~strcmp(filename, '')
        saveas(fig, filename);
        close(fig);
    end
    
    function doplot(startX, endX, showTitle)
        hold on;
%         for ii = 1 : statenumber
%             plot(X(startX : endX), allstates(startX:endX,ii,1), 'Color', PickColor(ii));
%         end
        for ii = 1 : lineNumber
            info = lines(ii);
            from = max(startX, info.StartX);
            to = min(endX, info.EndX);
            if from < to
                plot(X(from:to), allstates(from:to, info.StateId, 1), 'Color', info.Color, 'LineStyle', info.LineStyle);
            end
        end

        if showTitle
            title(strcat({'Lowest '},  num2str(statenumber), ' Energy States for M=', num2str(bits)));
        end
        
        ax1 = gca;
        xlabel('1/N');
        ylabel('E');
        if startX == 1
        else
            set(ax1, 'YAxisLocation', 'right');
        end
        
        set(ax1, 'XTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]);
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

