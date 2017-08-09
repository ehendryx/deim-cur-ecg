function [C,U,R,p,q] = CURfacQR(A,tol,truesvd)
% This function computes the CUR factorization of the matrix A
%
% Input:    A - an mxn matrix
%           tol - an error tolerance;
%                 If the SVD is found without incremental QR, this tolerance
%                 is used to determine the number of DEIM points to select.
%                 If incremental QR is used in approximating the SVD of A,
%                 this tolerance is divided by 10 prior to input into the
%                 incremental QR algorithm
%           truesvd (optional) - if truesvd = 1, then the SVD is computed directly from A;
%                                if truesvd = 0, then incremental QR is used;
%                                default = 0
%
% Output:   C - an mxk matrix such that C = A(:,q)
%           U - a kxk matrix such that U = C^I * A * R^I, where C^I
%               and R^I are the left and right inverses of C and R,
%               respectively
%           R - a kxn matrix such that R = A(p,:)
%           p - the DEIM-selected row indices
%           q - teh DEIM-selected column indices

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<3,
    truesvd = 0;
end

% Compute the SVD with or without incremental QR
if truesvd
    [V,S,W] = svd(full(A),'econ');
else
    tol1 = tol/10;
    [m,n] = size(A);
    
    % Apply incremental QR to A
    [V,R] = QRupdate_prealloc(full(A),tol1);
    
    % Approximate the SVD of A from the incremental QR output
    [V1,S,W] = svd(R,'econ');
    V = V*V1;    
end

% Isolate the sincgular values
s = diag(S);

% Determine rank of current SVD approximation
kk = size(V,2);

if truesvd
    % select the least index k  such that  s(k+1) < tol*s(1)
    j = 2;
    while (j <= kk) && (s(j) > tol*s(1))
        j = j+1;
    end
    
    if j > kk
        disp('failed to meet svd truncation tolerance  ')
    end
    
    k = j-1;
else   
    % Rely on incremental QR for approximate rank
    k = kk;
end

% Use rank-k SVD
V = V(:,1:k);
W = W(:,1:k);
S = S(1:k,1:k);

% Use DEIM to select columns and rows of A
[p] = Deim(V);
[q] = Deim(W);

% Form C, U, and R
C = A(:,q);
R = A(p,:);

[Q1,R1] = qr(full(C),0);
[Q2,R2] = qr(full(R'),0);

U = R1\(Q1'*A*Q2)/R2';




