function M = minors(X)
% minors
%   X   - Input matrix, square.
% 
% RETURN
%   M   - Output matrix, square.
%
% DESCRIPTION
%   Calculates the minors of the matrix. The minor M[i,j] is computed as
%   the determinant of the matrix X, whereas the i-th row and j-th column 
%   are deleted.
%
%   Copyright (C) 2013  Florian Raudies, 01/08/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

[yNum xNum] = size(X);
M = zeros(yNum,xNum);
for iy = 1:yNum,
    % Use symmetry.
    for ix = iy:xNum,
        C = X;
        C(iy,:) = [];
        C(:,ix) = [];
        d = det(C);
        M(iy,ix) = d;
        M(ix,iy) = d;
    end
end