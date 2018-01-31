function [landmark_set, landmarks, assign, dist] = ...
    select_landmark(data_struct, output_dir, num_landmark, num_core)
% 
% data_struct.name = 'csc_txt' | 'dense_mat' | 'sparse_mat'
% data_struct.files = {data_file, indices_file, indptr_file} | data_dir
% data_struct.data = sparse/dense

data_type = data_struct.name;
if strcmp(data_type, 'csc_txt')
    data_file = data_struct.files{1};
    indices_file = data_struct.files{2};
    indptr_file = data_struct.files{3};
elseif strcmp(data_type, 'dense_mat')
    [data_file, indices_file, indptr_file] = ...
        convert_dense_to_csc_matrix(data_struct.data, data_struct.files);
elseif strcmp(data_type, 'sparse_mat')
    [data_file, indices_file, indptr_file] = ...
        convert_sparse_to_csc_matrix(data_struct.data, data_struct.files);
end

exe_dir = data_struct.exe_dir;

command = sprintf('%s/landmark_selector -p %s -d %s -i %s -c %d -n %d -o %s',...
    exe_dir, indptr_file, data_file, indices_file, num_core, num_landmark, output_dir);
system(command);

% read landmarks
landmark_file = sprintf('%s/landmark.txt',output_dir);
[landmark_set, landmarks] = libsvmread(landmark_file);

% convert from c++ to matlab
landmark_set = landmark_set + 1;

assign_dist_file = sprintf('%s/assign_dist.txt', output_dir);
assign_dist = importdata(assign_dist_file);

assign = assign_dist(:,2);
% convert from c++ to matlab
assign = assign + 1;
dist = assign_dist(:,3);

function [data_file, indices_file, indptr_file] =...
    convert_dense_to_csc_matrix(X, data_dir)
% data: D x N matrix

if exist(data_dir,'dir') ~= 7
    mkdir(data_dir);
end

data_file = sprintf('%s/data.csv', data_dir);
indices_file = sprintf('%s/indices.csv', data_dir);
indptr_file = sprintf('%s/indptr.csv', data_dir);

[d, n] = size(X);
data = reshape(X,1,n*d);
indices=repmat(0:d-1, 1, n);
indptr=0:d:n*d;

dlmwrite(data_file, data,'delimiter', ',', 'precision',16);
dlmwrite(indices_file, indices,'delimiter', ',', 'precision',16);
dlmwrite(indptr_file, [d,n],'delimiter', ',', 'precision',16);
dlmwrite(indptr_file, indptr,'-append','delimiter', ',', 'precision',16);

function [data_file, indices_file, indptr_file] =...
    convert_sparse_to_csc_matrix(X, data_dir)
% data: D x N matrix
if exist(data_dir,'dir') ~= 7
    mkdir(data_dir);
end

data_file = sprintf('%s/data.csv', data_dir);
indices_file = sprintf('%s/indices.csv', data_dir);
indptr_file = sprintf('%s/indptr.csv', data_dir);

[d, n] = size(X);
data = [];
indices = [];
indptr = [0];
for i=1:n
   [rows, cols, vals] = find(X(:,i));
   data = [data, vals'];
   indices = [indices, rows'];
   count = indptr(end) + length(rows);
   indptr=[indptr, count];
end

dlmwrite(data_file, data,'delimiter', ',', 'precision',16);
dlmwrite(indices_file, indices, 'delimiter', ',', 'precision',16);
dlmwrite(indptr_file, [d,n],'delimiter', ',', 'precision',16);
dlmwrite(indptr_file, indptr,'-append', 'precision',16);