function p = Deim(U)
% This function carries out DEIM point selection using the columns of the matrix U.
%
% Input:    U - a full rank mxn matrix with m>=n
%
% Output:   p - an nx1 vector containing the indices of the DEIM-selected rows of U

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[m,n] = size(U);

p = zeros(n,1);

% Select the first DEIM point
[rho, p(1)] = max(abs(U(:,1)));

% Identify additional DEIM points, reading in additional columns of U one at a time
for j = 2:n,
    
    u = U(:,j);
    
    r = u - U(:,1:j-1)*(U(p(1:j-1),1:j-1)\u(p(1:j-1)));
    
    [rho,p(j)] = max(abs(r));   
end



