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
    
    % first go through midX to tot and connect lines.
    for curr = midX : tot
        n = 1/X(curr);
        states = lowest_energies(bits, n, statenumber);
        initDegenerate;
        
        if curr == midX
            allstates(midX) = states;
            highestE = states(statenumber, 1);
            lowestE = states(1, 1);
            initGap = highestE - lowestE;
        else
            build_dismat;
            res = find_match(dismat, statenumber);
            for j = 1 : statenumber
                allstates(curr, j) = states(res(j));
            end
            
            highestE = max(states(statenumber, 1), highestE);
            lowestE = min(states(1, 1), lowestE);
            currGap = highestE - lowestE;
            
            if currGap > 3 * initGap
                newPlotStart = curr;
            end
        end
        
        nextPos;
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
            allstates(curr, j) = states(res(j));
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
            f = allstates(p1, stateid);
        else
            f = (allstates(p2, stateid) - allstates(p1, stateid)) * (curr - p2) / (p2 - p1);
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
            for jj = 1 : statenumber
                if isDegenerate(jj)
                    dismat((ii-1)*statenumber + jj, 1) = norm(next(1:2) - states(jj,1:2));
                else
                    dismat((ii-1)*statenumber + jj, 1) = norm(next - states(jj));
                end
                
                dismat((ii-1)*statenumber + jj, 1) = ii;
                dismat((ii-1)*statenumber + jj, 1) = jj;
            end
        end
    end
    
end

function f=IsEqual(a, b)
    f = abs(a - b) < 1e-6;
end
        