function [k_nn,k_smallest_distances] = find_k_nearest_neighbors(to_match,library,k,distance)
% This function identifies the k nearest neighbors and corresponding k
% smallest neighbor distances
%
% Input:    to_match - a matrix of neighbors to search via k-nn with each
%                      row corresponding to a different point (neighbor)
%           library - a matrix containing points for which to_match
%                      neighbors are found; each row contains a point for
%                      k-nn analysis
%           k - number of nearest neighbors to find
%           distance - the distance measure of choice; the only current
%           option is distance = 'euclidean';
%
% Output:   k_nn - the indices of the k nearest neigbors listed in
%                  ascending order according to distance
%           k_smallest_distances - the corresponding distances to the
%                  neighors identified in k_nn


if strcmpi(distance,'euclidean')
    [k_nn,k_smallest_distances] = knnsearch(library,to_match,'K',k,'Distance',distance);
else % Could include other distances here in the future
    display('Distance measure is not an option. The only current option is ''euclidean''.')
    return
end


