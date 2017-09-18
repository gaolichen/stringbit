function [len, data] = norm_raw_data(bits, fNumber)
%norm_raw_data Reads raw data from norm data file.
%   Detailed explanation goes here
%   The norm matrix a symmetric matrix of dimension len * len
    file = strcat('norm', num2str(bits), '_', num2str(fNumber), '.dat');
    fid = fopen(file);
    line = fgetl(fid);
    cell = textscan(line, '%d');
%     display(cell)
    len = cell{1};
    
    data = zeros((len + 1) * len / 2, floor((bits + 1) / 2) + 1);
    count = 1;
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line)
            % do nothing.
        else
            cell = textscan(line, '%d');
            dim = size(cell{1});
%             display(cell{1})
%             display(dim)
            for i = 1 : dim(1)
                data(count, i) = cell{1}(i);
            end
            count = count + 1;
        end
    end
    
    fclose(fid);
end

