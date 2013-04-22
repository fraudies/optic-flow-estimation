function R = lsqConv3D(F, A, opt)
% lsqConv3D
%   F       - Input: 3D matrix with dimensions: yNum x xNum x tNum.
%   A       - Weight 1D kernel applied in all three dimension with kNum
%             entries.
%   opt     - Options with fields:
%
% RETURN
%   R       - components in following order:
%
%                 /r8 r5 r6 \      /r2\
%             A = |r5 r9 r7 |, b = |r3|, c=r1
%                 \r6 r7 r10/      \r4/
%
%             with
%
%             F ~= (x,y,t)^t A (x,y,t) + b^t (x,y,t) + c
%
%             The order of R corresponds to the base functions
%
%             B = {1, x, y, t, xy, xt, yt, x^2, y^2, t^2}
%
%             for each 2D position. R is a 4D matrix with dimensions:
%             yNum x xNum x tNum x 10.
%
%
%   Copyright (C) 2013  Florian Raudies, 01/04/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

[yNum xNum tNum] = size(F);
A       = A(:);
kNum    = length(A);

if ~isfield(opt,'SelT'), opt.SelT = [1 tNum-kNum+1]; end
SelT = opt.SelT;

% Check the size of filter kernel.
if ~mod(kNum, 2), error('Odd kernel size required!'); end

% Construct a basis for 3D coordinates.
kNum2 = (kNum-1)/2;
[Y X T] = ndgrid( -kNum2:+kNum2, -kNum2:+kNum2, -kNum2:+kNum2);
B = cat(4, ...
    ones(kNum,kNum,kNum), X, Y, T, X.*Y, X.*T, Y.*T, X.^2, Y.^2, T.^2);
bNum = size(B, 4); % Number of base vectors.

% Define the separable filter kernel in one dimension.
Ka = repmat(A(:)*A(:)',[1 1 kNum]).*repmat(reshape(A(:),[1 1 kNum]),...
    [kNum kNum 1]);

% Metric of base B(:,:,i) .* Ka .* B(:,:,j), i,j=1..kNum, and sum over i,j
M = sum(sum(sum(repmat(B,[1 1 1 1 bNum]) ...
             .* repmat(Ka,[1 1 1 bNum bNum]) ...
             .* permute(repmat(B,[1 1 1 1 bNum]),[1 2 3 5 4]))));
         
% Compute the inverse metric.
Minv = inv(squeeze(M));

% Build kernels for basis functions combined with weights a for the 
% neighborhood.
Xa = (-kNum2:+kNum2)';
K0 = A;
K1 = Xa.*A;
K2 = Xa.^2.*A;

% *************************************************************************
% Convolution of signal with basis functions
% *************************************************************************

% (i) Convolution in t dimension.
RconvT0 = convn(F, +reshape(K0,[1 1 kNum]), 'valid');
RconvT1 = convn(F, -reshape(K1,[1 1 kNum]), 'valid'); % '-' because of convn
RconvT2 = convn(F, +reshape(K2,[1 1 kNum]), 'valid');
RconvT0 = RconvT0(:,:,SelT(1):SelT(2));
RconvT1 = RconvT1(:,:,SelT(1):SelT(2));
RconvT2 = RconvT2(:,:,SelT(1):SelT(2));

% (ii) Convolution in x dimension.
RconvX0T0 = imfilter(RconvT0, K0');
RconvX1T0 = imfilter(RconvT0, K1');
RconvX2T0 = imfilter(RconvT0, K2');
RconvX0T1 = imfilter(RconvT1, K0');
RconvX1T1 = imfilter(RconvT1, K1');
RconvX0T2 = imfilter(RconvT2, K0');

% (iii) Convolution in y dimension.
Rconv = zeros(yNum,xNum,SelT(2)-SelT(1)+1,bNum);

Rconv(:,:,:,1)  = imfilter(RconvX0T0, K0);      % y^0 x^0 t^0

Rconv(:,:,:,2)  = imfilter(RconvX1T0, K0);      % y^0 x^1 t^0
Rconv(:,:,:,3)  = imfilter(RconvX0T0, K1);      % y^1 x^0 t^0
Rconv(:,:,:,4)  = imfilter(RconvX0T1, K0);      % y^0 x^0 t^1

Rconv(:,:,:,5)  = imfilter(RconvX1T0, K1)/2;    % y^1 x^1 t^0
Rconv(:,:,:,6)  = imfilter(RconvX1T1, K0)/2;    % y^0 x^1 t^1
Rconv(:,:,:,7)  = imfilter(RconvX0T1, K1)/2;    % y^1 x^0 t^1

Rconv(:,:,:,8)  = imfilter(RconvX2T0, K0);      % y^0 x^2 t^0
Rconv(:,:,:,9)  = imfilter(RconvX0T0, K2);      % y^2 x^0 t^0
Rconv(:,:,:,10) = imfilter(RconvX0T2, K0);      % y^0 x^0 t^2

% Apply the inverse metric
% M^-1 [(Q B) * F], where kNum2 has the form
% M^-1 = [m11 0   0   0   0   0   0   m18 m19 m110,
%         0   m22 0   0   0   0   0   0   0   0,
%         0   0   m33 0   0   0   0   0   0   0,
%         0   0   0   m44 0   0   0   0   0   0,
%         0   0   0   0   m55 0   0   0   0   0,
%         0   0   0   0   0   m66 0   0   0   0,
%         m18 0   0   0   0   0   m77 0   0   0,
%         m19 0   0   0   0   0   0   m88 0   0,
%         m10 0   0   0   0   0   0   0   0 m1010]

R = zeros(size(Rconv));
for i=1:bNum, R(:,:,:,i) = Rconv(:,:,:,i) * Minv(i,i); end

R(:,:,:,1)  = R(:,:,:,1)    + Rconv(:,:,:,8)*Minv(8,1) ...
                            + Rconv(:,:,:,9)*Minv(9,1) ...
                            + Rconv(:,:,:,10)*Minv(10,1);
R(:,:,:,8)  = R(:,:,:,8)    + Rconv(:,:,:,1)*Minv(8,1);     % X^2
R(:,:,:,9)  = R(:,:,:,9)    + Rconv(:,:,:,1)*Minv(9,1);     % Y^2
R(:,:,:,10) = R(:,:,:,10)   + Rconv(:,:,:,1)*Minv(10,1);    % T^2
