function f = orthnormal(A, positiveonly)
%     B = orth(A);
%     C = ctranspose(B);
%     D = C * A * ctranspose(C);
%     for i = 1 : size(C, 1)
%         if abs(D(i, i)) > 1e-6
%             C(i,:) = C(i,:) / sqrt(abs(D(i, i)));
%         end
%     end
%     
%     f = C;

    if nargin < 2
        % by default, dont return vec.
        positiveonly = false;
    end
    n = size(A, 1);
    S = zeros(n, n);
    f = zeros(n, n);
    nom=zeros(1, n);
    
    for i = 1 : n
        S(i, i) = 1;
        for j = 1 : i - 1
            if abs(nom(j)) > 1e-6
                a = S(i, :) * A * ctranspose(S(j, :));
                %if abs(a) > 1e-6
                    S(i, :) = S(i, :) - a * S(j, :) / nom(j);
                %end
            end
        end
        
        nom(i) = S(i, :) * A * ctranspose(S(i, :));
        if nom(i) > 1e-6
            f(i, :) = S(i, :) / sqrt(nom(i));
        elseif nom(i) <-1e-6 && ~positiveonly
            f(i, :) = S(i, :) / sqrt(-nom(i));
        else
            f(i, :) = S(i, :);
        end
    end
end