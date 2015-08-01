function [f, v]=lowest_energies(bits, n, statenumber, retVec)
    if nargin < 4
        % by default, dont return vec.
        retVec = false;
    end
    [Vec,D] = eigs(ham(bits, n), statenumber, 'sr');
    state_norm = real(diag(ctranspose(Vec)*metric(bits, n)*Vec));
    if retVec == false
        f = sortrows(horzcat(real(diag(D)),imag(diag(D)), state_norm), 1);
        v = {};
    else
        [f, indices] = sortrows(horzcat(real(diag(D)),imag(diag(D)), state_norm), 1);
        
        v = cell(1, statenumber);
        i = 1;
        w = size(Vec, 1);
        orders = 1:1:w;
        maxStates = 50;
        while i <= statenumber
            if i < statenumber && norm(f(i, 1) - f(i+1, 1), f(i, 2) - f(i+1, 2)) < 1e-6
                v{i} = {};
                v{i+1} = {}; 
                i = i + 2;
            else
                pos = indices(i);
                if (w > maxStates)
                    [~, stateIndices] = sortrows(abs(Vec(:,pos)), -1);
                    arr = Vec(stateIndices(1:maxStates),pos);
                    %v{i} = horzcat(arr, stateIndices(1:maxStates));
                    v{i} = sortrows(horzcat(arr, stateIndices(1:maxStates)), 2);
                else
                    v{i} = sortrows(horzcat(Vec(:,pos), orders'), 2);
                end
                i = i + 1;
            end
        end
    end
end