function normalized_data_matrix = data_matrix_beat_normalization(data_matrix)
% This function Z-normalizes each row of the given data matrix
% (resulting in a mean of 0 and standard deviation of 1)

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[m,n] = size(data_matrix);

% Subtract out mean from each row:
data_matrix = data_matrix - (1/n)*data_matrix*ones(n,1)*ones(1,n);

normalized_data_matrix = NaN(m,n);

% Divide each row by standard deviation:
std_dev = ones(1,m);
for i = 1:m
    std_dev(i) = std(data_matrix(i,:));
    normalized_data_matrix(i,:) = data_matrix(i,:)/std_dev(i);
end


