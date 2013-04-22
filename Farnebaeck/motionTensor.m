function T = motionTensor(R, gamma)
% motionTensor
%   R       - Minimum residual of the input signal projected to the basis
%             1,x,y,t,xy,xt,yt,x^2,y^2,t^2 which is a 4D matrix with the
%             dimensions: yNum x xNum x tNum x 10.
%   gamma   - Parameter that steers the influence between odd part and even 
%             part of the tensor component.
%
% RETURN
%   T       - Dimensions of the motion tensor: yNum x xNum x tNum x 3 x 3.
%
% DESCRIPTION
%   The motion tensor is defined by:
%
%       T = A * A^t + gamma b * b^t
%
%   with
%
%           /r8 r5 r6  \    /r2\
%       A = |r5 r9 r7  |, b=|r3|
%           \r6 r7 r10 /    \r4/
%
%   and r2,...,r10 are the residual components.
%
%   Copyright (C) 2013  Florian Raudies, 01/04/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

[yNum xNum tNum, ~] = size(R);
T = zeros(yNum,xNum,tNum,3,3);
T(:,:,:,1,1) = R(:,:,:,8).^2 + R(:,:,:,5).^2 + R(:,:,:,6).^2 ...
             + gamma * R(:,:,:,2).^2;
T(:,:,:,1,2) = R(:,:,:,5).*(R(:,:,:,8) + R(:,:,:,9)) ...
             + R(:,:,:,6).*R(:,:,:,7) ...
             + gamma * R(:,:,:,2).*R(:,:,:,3);
T(:,:,:,1,3) = R(:,:,:,6).*(R(:,:,:,8) + R(:,:,:,10)) ...
             + R(:,:,:,5).*R(:,:,:,7) ...
             + gamma * R(:,:,:,2).*R(:,:,:,4);
T(:,:,:,2,1) = T(:,:,:,1,2);
T(:,:,:,2,2) = R(:,:,:,5).^2 + R(:,:,:,9).^2 + R(:,:,:,7).^2 ...
             + gamma * R(:,:,:,3).^2;
T(:,:,:,2,3) = R(:,:,:,7).*(R(:,:,:,9) + R(:,:,:,10)) ...
             + R(:,:,:,5).*R(:,:,:,6) ...
             + gamma * R(:,:,:,3).*R(:,:,:,4);
T(:,:,:,3,1) = T(:,:,:,1,3);
T(:,:,:,3,2) = T(:,:,:,2,3);
T(:,:,:,3,3) = R(:,:,:,6).^2 + R(:,:,:,7).^2 + R(:,:,:,10).^2 ...
             + gamma * R(:,:,:,4).^2;
