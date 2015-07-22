function plot_states(bits, statenumber, points, ylims)
    if nargin < 3
        % by default, pick 100 points.
        points = 100;
    end
    if nargin < 4
        ylims = [];
    end
    
    minN = 1/1.5;
    X = 0:1/minN/points:1/minN;
    Y = zeros(length(X), statenumber);
    Yi = zeros(length(X), statenumber);
    Z = zeros(length(X), statenumber);
    slopeY = zeros(1, statenumber);
    slopeZ = zeros(1, statenumber);
    slopeYi = zeros(1, statenumber);
    complexLines = [];
    currComplexLines = [];
    inComplexLines = zeros(1, statenumber);
    
%     nextNorm = zeros(1, statenumber);
    
    for i = 1 : length(X)
        if X(i) == 0
            n = Inf;
        else
            n = 1/X(i);
        end
        states = lowest_energies(bits, n, statenumber);
        if i < 2
            for j = 1 : statenumber
                Y(i, j) = states(j, 1);
                Yi(i, j) = states(j, 2);
                Z(i, j) = states(j, 3);
            end
        else
            mat = zeros(statenumber * statenumber, 3);
            for j = 1 : statenumber
                for k = 1 : statenumber
                    mat((j-1)*statenumber + k, 1) = norm([states(k, 1) - slopeY(j), states(k, 2) - slopeYi(j),  states(k, 3) - slopeZ(j)]);
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
            end
        end
        
        for j = 1 : statenumber
            % if it's complex and not in inComplexLines 
            if NonZero(Yi(i, j)) && ~inComplexLines(j)
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
            
            if i == 1
                slopeY(j) = Y(i, j);
                slopeZ(j) = Z(i, j);
                slopeYi(j) = Yi(i, j);
%                 nextNorm(j) = 1.0;
%             elseif i == 2
            else
                slopeY(j) = 2 * Y(i, j) - Y(i - 1, j);
                slopeZ(j) = 2 * Z(i, j) - Z(i - 1, j);
                slopeYi(j) = 2 * Yi(i, j) - Yi(i - 1, j);
%                 nextNorm(j) = norm([slopeY(j) - Y(i, j), slopeZ(j) - Z(i, j), slopeYi(j) - Yi(i, j), X(i) - X(i-1)]);
%             else
%                 slopeY(j) = 3 * Y(i, j) - 3 * Y(i - 1, j) + Y(i - 2, j);
%                 slopeZ(j) = 3 * Z(i, j) - 3 * Z(i - 1, j) + Z(i - 2, j);
%                 slopeYi(j) = 3 * Yi(i, j) - 3 * Yi(i - 1, j) + Yi(i - 2, j);
%                 nextNorm(j) = norm([slopeY(j) - Y(i, j), slopeZ(j) - Z(i, j), slopeYi(j) - Yi(i, j), X(i) - X(i-1)]);
            end
        end
    end
    
    % move the remaining currComplines to complexLines
    for i = 1 : size(currComplexLines, 1)
        complexLines = cat(1, complexLines, currComplexLines(i,:));
    end
    
    % segments{i}: contains the segment should be excluded for state i.
    % segments{i}(2k-1): the beginning of the k-th exclude segment.
    % segments{i}(2k): the end of the k-th exclude segment
    segments = cell(1, statenumber);
    for i = 1 : statenumber
        segments{i} = 1:length(X);
    end
    for i = 1 : size(complexLines, 1)
        segments{complexLines(i, 2)} = cat(1, segments{complexLines(i, 2)}, [complexLines(i, 3) complexLines(i, 4)]);
    end
    % if no excluded segments, let the beginning of 
    for i = 1 : statenumber
        if size(segments{i}, 1) == 0
            segments{i} = length(X) + 1;
        end
    end
    
    D = {};
    figure(1); cla;
    hold on;
    for i = 1 : statenumber
        st = 1;
        j = 1;
        solidLine = IsPositive(Z(1, i));
        stop = 1;
        %for j = 2 : length(X)
        while j <= length(X)
            positive = IsPositive(Z(j, i));
            if j == segments{i}(stop) || solidLine ~= positive
                % plot the line
                if solidLine
                    plot(X(st:j), Y(st:j, i), 'Color', PickColor(i));
                else
                    plot(X(st:j), Y(st:j, i), 'Color', PickColor(i), 'LineStyle', ':');
                end
                
                if j == segments{i}(stop)
                    j = segments{i}(stop + 1);
                    st = j;
                    stop = stop + 2;
                    if j + 1 <= length(X)
                        solidLine = IsPositive(Z(j + 1, i));
                    end
                else
                    st = j;
                    solidLine = positive;
                end
            end
            j = j + 1;
        end
        
        if solidLine
            plot(X(st : length(X)), Y(st : length(X), i), 'Color', PickColor(i));
        else
            plot(X(st : length(X)), Y(st : length(X), i), 'Color', PickColor(i), 'LineStyle', ':');
        end
        D = cat(2, D, X, Z(:, i));
    end
    %figure;
    %plot(C{:});
    title(strcat({'Lowest '},  num2str(statenumber), ' Energy States for M=', num2str(bits)));
    ylabel('E');
    xlabel('1/N');
    ax1 = gca;
    set(ax1, 'XTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]);
    if length(ylims) == 4
        ylim([ylims(1) ylims(2)]);
    else
        ylim auto;
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

function f = NonZero(a)
    f = abs(a) > 1e-6;
end

    function f = IsZero(a)
        f = abs(a) < 1e-6;
    end

function f = IsPositive(a)
    f = a > 1e-6;
end

function c = PickColor(index)
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
        c = 'y';
    else
        c = 'w';
    end
end

    
  