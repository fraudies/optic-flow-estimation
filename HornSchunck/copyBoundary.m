function X = copyBoundary(X, len)
% copyBoundary
%   X       - nD Input matrix.
%   len     - Length of boundary which should be repeated.
%
% RETURN
%   X       - Output matrix (n-d), where len values are copied into the 
%             matrix 
%             * position (len+1) for the left boundary
%             * position (end-len) for the right boundary.
%             Unlike the function repeatBoundary here a boundary defined
%             inside the matrix len values from the existing boundary is
%             COPIED!
%
%   Copyright (C) 2013  Florian Raudies, 01/01/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

Dims    = size(X);
n       = ndims(X);
idx     = cell(1, n);
for id = 1:n,
    m = Dims(id);
    if m==1,
        idx{id} = 1;
    else
        idx{id} = [repmat(len+1,[1 len]) len+1:m-len repmat(m-len,[1 len])];
    end
end
X = X(idx{:});
