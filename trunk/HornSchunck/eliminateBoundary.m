function X = eliminateBoundary(X, len)
% eliminateBoundary
%   X       - nD matrix.
%   len     - 'len' is eliminated from both boundaries for all 
%             non-singleton dimensions.
%
% RETURN
%   X       - nD matrix. All non-singleton dimension are reduced by a
%             2*len in their length.
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
        idx{id} = (len+1):(m-len);
    end
end
X = X(idx{:});
