function [energy, norm, state] = get_largeN_ground_state(bit, minFermion)
%get_largeN_ground_stae Returns the ground state when N=Inf.
%   minFermion: the minimum ferminon number of the trace states to be
%   returned.

    [V,D] = eigs(ham(bit, Inf), 1, 'sr');
    w0 = V(:, end);
    energy = D(1, 1);
    
    vec = [];
    state=[];
    strTbl=importdata(strcat('state_structure_', num2str(bit), '.txt'));
    
    for i = 1 : size(w0, 1)
        if abs(w0(i)) > 1e-8 && strTbl.data(i, 2) >= minFermion
            vec = cat(1, vec, [i, w0(i), abs(w0(i))]);
            state=cat(1, state, [i, w0(i), abs(w0(i))]);
        else
            vec = cat(1, vec, [i, 0, 0]);
        end
    end
    
    amp=vec(1:end,2);
    norm=ctranspose(amp)*metric(bit, Inf)*amp;
end

