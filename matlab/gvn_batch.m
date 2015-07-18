
root = '~/gitroot/stringbit/matlab';
cd(root);
fid=fopen('ground_vanish_N.txt','w');

addpath('../..');
addpath('..');

folders ={'xi=0' 'xi=1' 'xi=5' 'xi=10' 'xi=n1' 'xi=n5' 'xi=n10'};
%folders ={'xi=0'};

cd('boson');
%fprintf(fid, 'size=%d %d\n', size(folders, 1), size(folder, 2));
fprintf('processing boson folder\n');
fprintf(fid, 'boson\n');
for i = 1 : size(folders, 2)
    cd(folders{i});
    fprintf('processing %s ...\n', folders{i});
    fprintf(fid, '%s\n', folders{i});
    fprintf('ground_vanish_N\n');
    fprintf(fid, 'ground_vanish_N\n');
    res = ground_vanish_N(11);
    for j = 1 : size(res, 1)
        fprintf(fid, '%2.0f %4.1f\n', res(j, 1), res(j, 2));
    end
    
    fprintf('ground_vanish_N2\n');
    fprintf(fid, 'ground_vanish_N2\n');
    res = ground_vanish_N2(11);
    for j = 1 : size(res, 1)
        fprintf(fid, '%2.0f %4.1f\n', res(j, 1), res(j, 2));
    end
    
    cd('..');
end

% cd('../fermion');
% fprintf('processing fermion folder\n');
% 
% for i = 1 : size(folders, 2)
%     cd(folders{i});
%     fprintf('processing %s ...\n', folders{i});
%     fprintf(fid, '%s\n', folders{i});
%     res = ground_vanish_N(11);
%     for j = 1 : size(res, 1)
%         fprintf(fid, '%2.0f %4.1f\n', res(j, 1), res(j, 2));
%     end
%     cd('..');
% end

fclose(fid);
