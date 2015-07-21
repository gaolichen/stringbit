function f = plot_states(bits, statenumber, points, ylims)
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
        
        for j = 1 : statenumber
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
    
    C = cell(1, statenumber);
    D = {};
    cnt = 1;
    for i = 1 : statenumber
        %C = cat(2, C, X, Y(:, i));
        st = 1;
        isComplex = Nonzero(Yi(1, i));
        for j = 2 : size(Yi, 1)
            isNonzero = Nonzero(Yi(j, i));
            if isComplex ~= isNonzero
                fprintf('%d\n', j);
                C{1, cnt} = X(st : j - 1);
                cnt = cnt + 1;
                C{1, cnt} = Y(st : j - 1, i)';
                cnt = cnt + 1;
                if isComplex
                    %C = cat(2, C, X(st : j - 1), Y(st : j - 1, i)', '--');
                    C{1, cnt} = '--';
                    cnt = cnt + 1;
                else
                    %C = cat(2, C, X(st : j - 1), Y(st : j - 1, i)');
                end
                isComplex = isNonzero;
                st = j;
            end
        end
        if isComplex
            C = cat(2, C, X(st : size(X, 1)), Y(st : size(Y, 1), i)', '--');
        else
            C = cat(2, C, X(st : size(X, 1)), Y(st : size(Y, 1), i)');
        end
        D = cat(2, D, X, Z(:, i));
    end
    f = C;
%     figure;
%     plot(C{:});
%     title(strcat({'Lowest '},  num2str(statenumber), ' Energy States for M=', num2str(bits)));
%     ylabel('E');
%     xlabel('1/N');
%     ax1 = gca;
%     set(ax1, 'XTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]);
%     if length(ylims) == 4
%         ylim([ylims(1) ylims(2)]);
%     else
%         ylim auto;
%     end
%     
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

function f = Nonzero(a)
    f = abs(a) > 1e-6;
end

    
  