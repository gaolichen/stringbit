function plot_states2(bits, statenumber, points, ylims)
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
    
    for i = 1 : length(X)
        if X(i) == 0
            n = Inf;
        else
            n = 1/X(i);
        end
        states = lowest_energies2(bits, n, statenumber);
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
%                     mat((j-1)*statenumber + k, 1) = norm([states(k, 1) - slopeY(j), states(k, 2) - slopeYi(j)]);
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
        
        for j = 1 : statenumber
            if i == 1
                slopeY(j) = Y(i, j);
                slopeZ(j) = Z(i, j);
                slopeYi(j) = Yi(i, j);
%             elseif i < 10
            else
                slopeY(j) = 2 * Y(i, j) - Y(i - 1, j);
                slopeZ(j) = 2 * Z(i, j) - Z(i - 1, j);
                slopeYi(j) = 2 * Yi(i, j) - Yi(i - 1, j);
%             else
%                 slopeY(j) = 3 * Y(i, j) - 3 * Y(i - 1, j) + Y(i - 2, j);
%                 slopeZ(j) = 3 * Z(i, j) - 3 * Z(i - 1, j) + Z(i - 2, j);
%                 slopeYi(j) = 3 * Yi(i, j) - 3 * Yi(i - 1, j) + Yi(i - 2, j);
            end
        end
        
%         if i > 0 && i < 10
%              fprintf('=============== i = %d\n', i);
% %             for j = 1 : size(mat, 1)
% %                 fprintf('%f %d %d\n', mat(j, 1), int32(mat(j, 2)), int32(mat(j, 3)));
% %             end
% %             for j = 1 : statenumber
% %                 fprintf('%d ', res(j));
% %             end
% %             fprintf('\n');
%             for j = 1 : statenumber
%                 fprintf('%f %f %f\n', Y(i, j), Yi(i, j), Z(i, j));
%             end
%         end
    end
    
    C = {};
    D = {};
    for i = 1 : statenumber
        C = cat(2, C, X, Y(:, i));
        D = cat(2, D, X, Z(:, i));
    end
    figure;
    plot(C{:});
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
    
    figure;
    plot(D{:});
    title(strcat({'Norm of Lowest '}, num2str(statenumber), ' Energy States for M=', num2str(bits)));
    ylabel('Norm');
    xlabel('1/N');
    ax2 = gca;
    set(ax2, 'XTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]);
    if length(ylims) == 4
        ylim([ylims(3) ylims(4)]);
    else
        ylim auto;
    end
end