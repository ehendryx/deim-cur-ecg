function full_class = synthetic_class_counter(q)
% This function counts the number of representatives selected via DEIM CUR
% for the 12 synthetic classes
% 
% Input:    q - a vector of DEIM-selected column indices
%
% Output:   full_class - a 12x1 vector with each entry containing the
%                       number of identified representatives for the 
%                       corresponding class

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

full_class = zeros(12,1);
for t = 1:length(q)
    if q(t) <= 500
        full_class(1) = full_class(1) + 1;
    elseif 500 < q(t) && q(t) <= 1000
        full_class(2) = full_class(2) + 1;
    elseif 1000 < q(t) && q(t) <= 1500
        full_class(3) = full_class(3) + 1;
    elseif 1500 < q(t) && q(t) <= 2000
        full_class(4) = full_class(4) + 1;
    elseif 2000 < q(t) && q(t) <= 2500
        full_class(5) = full_class(5) + 1;
    elseif 2500 < q(t) && q(t) <= 3000
        full_class(6) = full_class(6) + 1;
    elseif 3000 < q(t) && q(t) <= 3500
        full_class(7) = full_class(7) + 1;
    elseif 3500 < q(t) && q(t) <= 4000
        full_class(8) = full_class(8) + 1;
    elseif 4000 < q(t) && q(t) <= 4500
        full_class(9) = full_class(9) + 1;
    elseif 4500 < q(t) && q(t) <= 5000
        full_class(10) = full_class(10) + 1;
    elseif 5000 < q(t) && q(t) <= 5500
        full_class(11) = full_class(11) + 1;
    elseif 5500 < q(t) && q(t) <= 6000
        full_class(12) = full_class(12) + 1;
    end
end