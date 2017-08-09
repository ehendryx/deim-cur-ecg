function [V,R] = QRupdate_prealloc(A,tol)
% This function applies incremental QR to the matrix A
%
% Input:    A - an mxn matrix
%           tol - the threshold used in the deflation step of the algorithm
%
% Output:   V - an mxk orthogonal matrix, and
%           R - a kxn matrix such that V*R is a rank-k QR factorization of A

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[m,n] = size(A);

% We choose our initial k to be 1, but this can take on different values if desired, e.g. k = min(20,n);
k = min(1,n);

% Initialize QR factorization
f = A(:,1);
rho = norm(f);  v = f/rho;
V = zeros(m,k);
R = zeros(k,n);
V(:,1) = v;  R(1,1) = rho;
s = zeros(n,1);

%    % Initial QR factorization to be included for starting k > 1
%    for j = 2:k,
%
%        a = A(:,j);
%        r = V(:,1:j-1)'*a;
%        f = a - V(:,1:j-1)*r;
%            % DGKS reorthogonalization
%            c = V(:,1:j-1)'*f;
%            f = f - V(:,1:j-1)*c;
%            r = r + c;
%        rho = norm(f);
%        v = f/rho;
%        V(:,j) = v;
%        R(1:j,j) = [r; rho];
%
%    end


% Initialize the squared row norm (held in s)
for i = 1:k,
    s(i) = (norm(R(i,:)))^2;
end


% Continue to form the QR factorization, reading in one column at a time 
% and deflating when a column contributes little information

j = k+1;
while(j <= n)
    
    a = A(:,j);
    r = V(:,1:k)'*a;
    f = a - V(:,1:k)*r;
        % DGKS reorthogonalization
        c = V(:,1:k)'*f;
        f = f - V(:,1:k)*c;
        r = r + c;
    rho = norm(f);
    v = f/rho;
    V(:,k+1) = v;
    R(1:k+1,j) = [r; rho];
    
    % Update the squared row norms and squared Frobenius norm (in FnormR)
    s(k+1) = rho^2;
    for i = 1:k,
        s(i) = s(i) + r(i)^2;
    end
    FnormR = sum(s);
    
    [sigma, i_min] = min(s(1:k+1));
    
    if (sigma > (tol^2)*(FnormR - s(i_min))),
        % No deflation
        k = k + 1;
    else
        % Use deflation
        if(i_min < k+1),
            % Row interchange R and column interchange V
            R(i_min,:) = R(k+1,:);
            V(:,i_min) = V(:, k+1);
            s(i_min) = s(k+1);
        end
        % Delete last row of R and last column of V
        V(:,k+1) = zeros(m,1);
        R(k+1,:) = zeros(1,n);
    end
    j = j+1;
    
    
end
V = V(:,1:k);
R = R(1:k,:);


