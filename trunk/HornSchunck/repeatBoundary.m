function X = repeatBoundary(X, len)
% repeatBoundary
%   X       - nD matrix.
%   len     - Length of boundary which should be repeated.
%
% RETURN
%   X       - nD matrix, whereas all non-singelton dimensions are enlarged 
%             by 2*len.
%
%   Copyright (C) 2013  Florian Raudies, 01/02/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

Dims    = size(X);
n       = ndims(X);
idx     = cell(1, n);
for id = 1:n,
    m = Dims(id);
    if m==1,
        idx{id} = 1;
    else
        idx{id} = [ones([1 len]) 1:m repmat(m,[1 len])];
    end
end
X = X(idx{:});
