function plot_states(bits, statenumber, points, maxX, filename)
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
    xInterval = round(0.02/(maxX) * points);
    X = 0:1/minN/points:1/minN;
    xStart = round((1/(1.5*bits)/maxX) * points);
    Y = zeros(length(X), statenumber);
    Yi = zeros(length(X), statenumber);
    Z = zeros(length(X), statenumber);
    slopeY = zeros(1, statenumber);
    slopeYi = zeros(1, statenumber);
    slopeZ = zeros(1, statenumber);
    lastReal = ones(1, statenumber);
    complexLines = [];
    currComplexLines = [];
    inComplexLines = zeros(1, statenumber);

    for j = 1 : statenumber
        lastReal(j) = xStart;
    end
    
    for i = xStart : length(X)
        if X(i) == 0
            n = Inf;
        else
            n = 1/X(i);
        end
        states = lowest_energies(bits, n, statenumber);
        if i < xStart + 1
            for j = 1 : statenumber
                Y(i, j) = states(j, 1);
                Yi(i, j) = states(j, 2);
                Z(i, j) = states(j, 3);
            end
        else
            mat = zeros(statenumber * statenumber, 3);
            for j = 1 : statenumber
                for k = 1 : statenumber
                    mat((j-1)*statenumber + k, 1) = norm([states(k, 1) - slopeY(j), states(k, 2) - slopeYi(j), states(k, 3) - slopeZ(j)]);
                    mat((j-1)*statenumber + k, 2) = j;
                    mat((j-1)*statenumber + k, 3) = k;
                end
            end
            res = find_match(mat, statenumber);
            for j = 1 : statenumber
                Y(i, j) = states(res(j), 1);
                Yi(i, j) = states(res(j), 2);
                Z(i, j) = states(res(j), 3);
            end
        end
        
        % check the complexLines.
        for j = size(currComplexLines, 1) : -1 : 1
            currComplexLines(j, 4) = i;
            line1 = currComplexLines(j, 1);
            line2 = currComplexLines(j, 2);
            if IsZero(Y(i, line1) - Y(i, line2)) && IsZero(Yi(i, line1) + Yi(i, line2))
                % do nothing
            else
                % the current complexLine ends, move it to complexLines.
                complexLines = cat(1, complexLines, currComplexLines(j,:));
                currComplexLines(j,:) = [];
                inComplexLines(line1) = 0;
                inComplexLines(line2) = 0;
                lastReal(line1) = i;
                lastReal(line2) = i;
            end
        end
        
        for j = 1 : statenumber
            % if it's complex and not in inComplexLines 
            if ~IsZero(Yi(i, j)) && ~inComplexLines(j)
                for k = j + 1 : statenumber
                    % If they conjugate to each other.
                    if IsZero(Yi(i, k) + Yi(i, j)) && IsZero(Y(i, k) - Y(i, j))
                        if inComplexLines(k)
                            fprintf('dotted line conflicts...');
                        end
                        inComplexLines(k) = 1;
                        inComplexLines(j) = 1;
                        currComplexLines = cat(1, currComplexLines, [j k i i]);
                        break;
                    end
                end
            end
            
            if inComplexLines(j)
                slopeY(j) = Y(i, j);
                slopeYi(j) = Yi(i, j);
                slopeZ(j) = Z(i, j);
            elseif i - lastReal(j) + 1 < 2 * xInterval
                slopeY(j) = Y(i, j);
                slopeYi(j) = Yi(i, j);
                slopeZ(j) = Z(i, j);
            else
                slopeY(j) = 2 * Y(i + 1 - xInterval, j) - Y(i + 1 - 2 * xInterval, j);
                slopeYi(j) = 2 * Yi(i + 1 - xInterval, j) - Yi(i + 1 - 2 * xInterval, j);
                slopeZ(j) = 2 * Z(i + 1 - xInterval, j) - Z(i + 1 - 2 * xInterval, j);
                if norm(slopeY(j) - Y(i, j), slopeYi(j) - Yi(i, j)) > 90 * X(2)
                    slopeY(j) = Y(i, j);
                    slopeYi(j) = Yi(i, j);
                    slopeZ(j) = Z(i, j);
                end
            end
        end
    end
    
    % process the beginning part. 
    for i = xStart - 1 : -1 : 1
        if X(i) == 0
            n = Inf;
        else
            n = 1/X(i);
        end
        states = lowest_energies(bits, n, statenumber);
        for j = 1 : statenumber
            for k = 1 : statenumber
                mat((j-1)*statenumber + k, 1) = norm([states(k, 1) - Y(i + 1, j), states(k, 2) - Yi(i + 1, j)]);
                mat((j-1)*statenumber + k, 2) = j;
                mat((j-1)*statenumber + k, 3) = k;
            end
        end
        
        res = find_match(mat, statenumber);
        for j = 1 : statenumber
            Y(i, j) = states(res(j), 1);
            Yi(i, j) = states(res(j), 2);
            Z(i, j) = states(res(j), 3);
        end
    end
    
    % move the remaining currComplines to complexLines
    for i = 1 : size(currComplexLines, 1)
        complexLines = cat(1, complexLines, currComplexLines(i,:));
    end
    
%     for i = 1 : size(complexLines, 1)
%         fprintf('%d %d %f %f\n', complexLines(i, 1), complexLines(i, 2), X(complexLines(i, 3)), X(complexLines(i, 4)));
%     end
    
    % segments{i}: contains the segment should be excluded for state i.
    % segments{i}(2k-1): the beginning of the k-th exclude segment.
    % segments{i}(2k): the end of the k-th exclude segment
    segments = cell(1, statenumber);
    for i = 1 : statenumber
        segments{i} = [];
    end
    for i = 1 : size(complexLines, 1)
        segments{complexLines(i, 2)} = [segments{complexLines(i, 2)}, [complexLines(i, 3) complexLines(i, 4)]];
    end
    % Append length(X) + 1 to the end for every line.
    for i = 1 : statenumber
        segments{i} = [segments{i}, length(X) + 1];
    end
    
%     for i = 1 : statenumber
%         fprintf('%d: ', i);
%         for j = 1 : size(segments{i}, 2)
%             if segments{i}(j) <= length(X)
%                 fprintf('%f ', X(segments{i}(j)));
%             else
%                 fprintf('end ');
%             end
%         end
%         fprintf('\n');
%     end
    
%    D = {};
    fig = figure; cla;
    hold on;
    for i = 1 : statenumber
        st = 1;
        j = 1;
        
        % lineType:
        % 0 for unset; 
        % 1 for solid line; 
        % 2 for dotted line
        lineType = 0;
        k = 1;
        while j <= length(X)
            if j == segments{i}(2 * k - 1);
                % when j is the start point of a excluded segments
                % plot here.
                if lineType ~= 0
                    plot(X(st:j), Y(st:j, i), 'Color', PickColor(i), 'LineStyle', PickLineStyle(lineType, X(j)), 'LineWidth', PickLineWidth(lineType));
                end
                
                j = segments{i}(2 * k);
                k = k + 1;
                % reset lineType
                lineType = 0;
            elseif lineType == 0
                % if lineType is not set.
                st = j;
                lineType = GetLineType(Z(j, i));
            elseif lineType ~= GetLineType(Z(j, i))
                % if the current line type does not match.
                % plot
                plot(X(st:j), Y(st:j, i), 'Color', PickColor(i), 'LineStyle', PickLineStyle(lineType, X(j)), 'LineWidth', PickLineWidth(lineType));
                lineType = GetLineType(Z(j, i));
                st = j;
            end
            j = j + 1;
        end
        
        % plot the last part.
        if lineType ~= 0 && st <= length(X)
            plot(X(st:length(X)), Y(st:length(X), i), 'Color', PickColor(i), 'LineStyle', PickLineStyle(lineType, 0),'LineWidth', PickLineWidth(lineType));
        end
%        D = cat(2, D, X, Z(:, i));
    end
    title(strcat({'Lowest '},  num2str(statenumber), ' Energy States for M=', num2str(bits)));
    ylabel('E');
    xlabel('1/N');
    ax1 = gca;
    set(ax1, 'XTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]);
    xlim([0, maxX]);
    ylim auto;
    
    if ~strcmp(filename, '')
        saveas(fig, filename);
        close(fig);
    end
    
%     figure;
%     plot(D{:});
%     title(strcat({'Norm of Lowest '}, num2str(statenumber), ' Energy States for M=', num2str(bits)));
%     ylabel('Norm');
%     xlabel('1/N');
%     ax2 = gca;
%     set(ax2, 'XTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]);
%     if length(ylims) == 4
%         ylim([ylims(3) ylims(4)]);
%     else
%         ylim auto;
%     end    
end

function f = IsZero(a)
    f = abs(a) < 1e-6;
end

function f = GetLineType(norm)
    % for positive norm, return type 1 (solid line)
    % for non-positive norm, return type 2 (dotted line)
    if norm > 1e-6
        f = 1;
    else
        f = 2;
    end
end

function f = PickLineStyle(lineType, x)
    if lineType == 1
        f = '-';
    elseif lineType == 2
        f = ':';
    else
        fprintf('Unexpected line type: %d x=%f\n', lineType, x);
        f = ':';
    end
end

function f = PickLineWidth(lineType)
    if lineType == 1
        f = 0.5;
    elseif lineType == 2
        f = 1;
    else
        fprintf('Unexpected line type: %d\n', lineType);
        f = 0.5;
    end
end

function c = PickColor(n)
    index = mod(n, 10);
    if index == 0
        index = 10;
    end
    if index == 1
        c = 'r';
    elseif index == 2
        c = 'k';
    elseif index == 3
        c = 'g';
    elseif index == 4
        c = 'b';
    elseif index == 5
        c = 'm';
    elseif index == 6
        c = 'c';
    elseif index == 7
        c = [0.7 0.7 0.7];
    elseif index == 8
        c = [0.7 0.4 1];
    elseif index == 9
        c = [0.6 0 0.3];
    else
        c = [1 0.5 0];
    end
end

    
  