function [Dx Dy] = estimateOpticFlow2D(ImgSeq, opt)
% estimateOpticFlow2D
%   ImgSeq  - Image sequence as a cube with dimensions: 
%             height x width x frames.
%   opt     - Struture with options:
%             * W     - Weight matrix for that sums constraints in the 
%                       local neighborhood.
%             * DiffX - Kernel that approximtes the computation of the
%                       partial derivative in x.
%             * DiffY - Ditto for y.
%             * DiffT - Ditto for t.
%
% RETURN
%   Dx      - X-component of flow and its spatio-temporal derivatives.
%   Dy      - Y-component of flow and its spatio-temporal derivatievs.
%             Computes flow for the center frame of the image sequence,
%             thus, Dx and Dy have the dimensions: height x width x 4.
%             The four stands for the two flow components and its three
%             partial derivatives (in x, y, and t), respectively.
%
% DESCRIPTION
%   An implementation of (rigorous condition) 
%   Nagel, H.-H. (1987). On the estimation of optical flow: relations 
%       between different approaches and some new results. Artificial 
%       Intelligence 33, 299-324.
%
%   Copyright (C) 2013  Florian Raudies, 01/06/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for paraemters of the method.
if nargin<2,                opt         = struct();                     end
if ~isfield(opt,'W'),       opt.W       = fspecial('gaussian',[9 9],2); end
if ~isfield(opt,'DiffX'),   opt.DiffX   = [-1 8 0 -8 1]/12;             end
if ~isfield(opt,'DiffY'),   opt.DiffY   = opt.DiffX';                   end
if ~isfield(opt,'DiffT'),   opt.DiffT   = reshape([-1 0 1]/2,[1 1 3]);  end
% Retrieve default parmeters.
W       = opt.W;
DiffX   = opt.DiffX;
DiffY   = opt.DiffY;
DiffT   = opt.DiffT;
ktNum   = size(DiffT,3);
yNum    = size(ImgSeq,1);
xNum    = size(ImgSeq,2);
tNum    = size(ImgSeq,3);
% Check if the provided sequence contains enough frames.
if tNum<(6*(ktNum-1)/2+1), 
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], 6*(ktNum-1)/2+1, tNum);
end
% Select the center number of 6*(ktNum-1)/2+1 frames.
center  = ceil(tNum/2);
SelT    = center + ( -3*(ktNum-1)/2 : +3*(ktNum-1)/2 );
ImgSeq  = ImgSeq(:,:,SelT);
% Compute 1st-order derivatives for the temporally valid part.
tNum    = size(ImgSeq,3);
ValidT  = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
F.x     = imfilter(ImgSeq(:,:,ValidT), DiffX, 'same', 'replicate');
F.y     = imfilter(ImgSeq(:,:,ValidT), DiffY, 'same', 'replicate');
F.t     = convn(ImgSeq, DiffT, 'valid');
% Compute 2nd-order derivatives for the temporally valid part.
tNum    = size(F.x,3);
ValidT  = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
F.xx    = imfilter(F.x(:,:,ValidT), DiffX, 'same', 'replicate');
F.xy    = imfilter(F.x(:,:,ValidT), DiffY, 'same', 'replicate');
F.yy    = imfilter(F.y(:,:,ValidT), DiffY, 'same', 'replicate');
F.xt    = convn(F.x, DiffT, 'valid');
F.yt    = convn(F.y, DiffT, 'valid');
F.tt    = convn(F.t, DiffT, 'valid');
% Compute 3rd-order derivatives for the temporally valid part.
tNum    = size(F.xx,3);
ValidT  = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
F.xxx   = imfilter(F.xx(:,:,ValidT), DiffX, 'same', 'replicate');
F.xxy   = imfilter(F.xx(:,:,ValidT), DiffY, 'same', 'replicate');
F.xxt   = convn(F.xx, DiffT, 'valid');
F.xyy   = imfilter(F.yy(:,:,ValidT), DiffX, 'same', 'replicate');
F.yyy   = imfilter(F.yy(:,:,ValidT), DiffY, 'same', 'replicate');
F.yyt   = convn(F.yy, DiffT, 'valid');
F.xtt   = imfilter(F.tt(:,:,ValidT), DiffX, 'same', 'replicate');
F.ytt   = imfilter(F.tt(:,:,ValidT), DiffY, 'same', 'replicate');
F.ttt   = convn(F.tt, DiffT, 'valid');
F.xyt   = convn(F.xy, DiffT, 'valid');
% Select the center frame for 1st and 2nd order derivatives.
center1 = ceil(size(F.x,3)/2);
F.x = squeeze(F.x(:,:,center1));
F.y = squeeze(F.y(:,:,center1));
F.t = squeeze(F.t(:,:,center1));
center2 = ceil(size(F.xx,3)/2);
F.xx = squeeze(F.xx(:,:,center2));
F.xy = squeeze(F.xy(:,:,center2));
F.yy = squeeze(F.yy(:,:,center2));
F.xt = squeeze(F.xt(:,:,center2));
F.yt = squeeze(F.yt(:,:,center2));
F.tt = squeeze(F.tt(:,:,center2));
% Apply smoothing to the computed derivatives. Note the original method did
% describe such smoothing.
Name = fieldnames(F);
nNum = length(Name);
for iName = 1:nNum,
    F.(Name{iName}) = imfilter(F.(Name{iName}), W, 'same', 'replicate');
end
% Allocate memory for the two flow components and their derivatives for the
% center frame.
Dx  = zeros(yNum, xNum, 4);
Dy  = zeros(yNum, xNum, 4);
for ix = 1:xNum,
    for iy = 1:yNum,
        % The rigorous condition provides 20 constraint equations. The
        % right-hand side of these equations is stored in B and the 
        % left-hand-side in the matrix A.
        B = [F.t(iy,ix) F.xt(iy,ix) F.yt(iy,ix) F.tt(iy,ix) 1/2*F.xxt(iy,ix) ...
             F.xyt(iy,ix) F.xtt(iy,ix) 1/2*F.yyt(iy,ix) F.ytt(iy,ix) ...
             1/2*F.ttt(iy,ix) 0 0 0 0 0  0 0 0 0 0]';
        A = [F.x(iy,ix)      0               0               0               F.y(iy,ix)      0               0               0             ;...
             F.xx(iy,ix)     F.x(iy,ix)      0               0               F.xy(iy,ix)     F.y(iy,ix)      0               0             ;...
             F.xy(iy,ix)     0               F.x(iy,ix)      0               F.yy(iy,ix)     0               F.y(iy,ix)      0             ;...
             F.xt(iy,ix)     0               0               F.x(iy,ix)      F.yt(iy,ix)     0               0               F.y(iy,ix)    ;...
             F.xxx(iy,ix)/2  F.xx(iy,ix)     0               0               F.xxy(iy,ix)/2  F.xy(iy,ix)     0               0             ;...
             F.xxy(iy,ix)    F.xy(iy,ix)     F.xx(iy,ix)     0               F.xyy(iy,ix)    F.yy(iy,ix)     F.xy(iy,ix)     0             ;...
             F.xxt(iy,ix)    F.xt(iy,ix)     0               F.xx(iy,ix)     F.xyt(iy,ix)    F.yt(iy,ix)     0               F.xy(iy,ix)   ;...             
             F.xyy(iy,ix)/2  0               F.xy(iy,ix)     0               F.yyy(iy,ix)/2  0               F.yy(iy,ix)     0             ;...
             F.xyt(iy,ix)    0               F.xt(iy,ix)     F.xy(iy,ix)     F.yyt(iy,ix)    0               F.yt(iy,ix)     F.yy(iy,ix)   ;...
             F.xtt(iy,ix)/2  0               0               F.xt(iy,ix)     F.ytt(iy,ix)/2  0               0               F.yt(iy,ix)   ;...
             0               F.xxx(iy,ix)/2  0               0               0               F.xxy(iy,ix)/2  0               0             ;...
             0               F.xxy(iy,ix)    F.xxx(iy,ix)/2  0               0               F.xyy(iy,ix)    F.xxy(iy,ix)/2  0             ;...
             0               F.xxt(iy,ix)    0               F.xxx(iy,ix)/2  0               F.xyt(iy,ix)    0               F.xxy(iy,ix)/2;...
             0               F.xyy(iy,ix)/2  F.xxy(iy,ix)    0               0               F.yyy(iy,ix)/2  F.xyy(iy,ix)    0             ;...
             0               F.xyt(iy,ix)    F.xxt(iy,ix)    F.xxy(iy,ix)    0               F.yyt(iy,ix)    F.xyt(iy,ix)    F.xyy(iy,ix)  ;...
             0               F.xtt(iy,ix)/2  0               F.xxt(iy,ix)    0               F.ytt(iy,ix)/2  0               F.xyt(iy,ix)  ;...
             0               0               F.xyy(iy,ix)/2  0               0               0               F.yyy(iy,ix)/2  0             ;...
             0               0               F.xyt(iy,ix)    F.xyy(iy,ix)/2  0               0               F.yyt(iy,ix)    F.yyy(iy,ix)/2;...
             0               0               F.xtt(iy,ix)/2  F.xyt(iy,ix)    0               0               F.ytt(iy,ix)/2  F.yyt(iy,ix)  ;...
             0               0               0               F.xtt(iy,ix)/2  0               0               0               F.ytt(iy,ix)/2];        
        % Compute the least squares solution for the flow and its
        % derivatives.
        X = - pinv(A) * B;
        Dx(iy,ix,:) = X(1:4);
        Dy(iy,ix,:) = X(5:8);
    end
end
