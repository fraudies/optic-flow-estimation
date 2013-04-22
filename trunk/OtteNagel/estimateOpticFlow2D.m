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
%   An implementation of the idea of 
%   Otte, M. and Nagel, H.-H. (1994). Optical flow estimation: advances and 
%       comparisons. Jan-Olof Eklundh (Ed.) ECCV'94, LNCS 800, 51-60.
%
%   Copyright (C) 2013  Florian Raudies, 01/06/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for paraemters of the method.
if nargin<2,                opt         = struct(); end
if ~isfield(opt,'W'),       opt.W       = fspecial('gaussian',[9 9],2); end
if ~isfield(opt,'DiffX'),   opt.DiffX   = [-1 8 0 -8 1]/12; end
if ~isfield(opt,'DiffY'),   opt.DiffY   = opt.DiffX'; end
if ~isfield(opt,'DiffT'),   opt.DiffT   = reshape([-1 0 1]/2,[1 1 3]); end
% Retrieve default parmeters to be saved in their own variables.
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
% Select the central number of 6*(ktNum-1)/2+1 frames.
center  = ceil(tNum/2);
SelT    = center + ( -3*(ktNum-1)/2 : +3*(ktNum-1)/2 );
ImgSeq  = ImgSeq(:,:,SelT);
% Compute 1st-order the derivatives for the temporally valid part.
tNum    = size(ImgSeq,3);
ValidT  = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
F.x     = imfilter(ImgSeq(:,:,ValidT), DiffX, 'same', 'replicate');
F.y     = imfilter(ImgSeq(:,:,ValidT), DiffY, 'same', 'replicate');
F.t     = convn(ImgSeq, DiffT, 'valid');
% Compute 2nd-order the derivatives for the temporally valid part.
tNum    = size(F.x,3);
ValidT  = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
F.xx    = imfilter(F.x(:,:,ValidT), DiffX, 'same', 'replicate');
F.xy    = imfilter(F.x(:,:,ValidT), DiffY, 'same', 'replicate');
F.yy    = imfilter(F.y(:,:,ValidT), DiffY, 'same', 'replicate');
F.xt    = convn(F.x, DiffT, 'valid');
F.yt    = convn(F.y, DiffT, 'valid');
F.tt    = convn(F.t, DiffT, 'valid');
% Compute 3rd-order the derivatives for the temporally valid part.
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
% Select the center frame for the 1st and 2nd order derivatives.
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
% Apply smoothing to the computed derivatives -- this was not described in 
% the original method.
Name = fieldnames(F);
nNum = length(Name);
for iName = 1:nNum,
    F.(Name{iName}) = imfilter(F.(Name{iName}), W, 'same', 'replicate');
end
% Allocate memory for the computed flow and its derivatives for the center 
% frame.
Dx  = zeros(yNum, xNum, 4);
Dy  = zeros(yNum, xNum, 4);
for ix = 1:xNum,
    for iy = 1:yNum,
        % Load data into auxiliary variables fx, fy, ft, fxx, fyy, ...
        % 1st order
        fx = F.x(iy,ix);
        fy = F.y(iy,ix);
        ft = F.t(iy,ix);
        % 2nd order
        fxx = F.xx(iy,ix);
        fyy = F.yy(iy,ix);
        ftt = F.tt(iy,ix);
        fxt = F.xt(iy,ix);
        fyt = F.yt(iy,ix);
        fxy = F.xy(iy,ix);
        % 3rd order
        fxxx = F.xxx(iy,ix);
        fxxy = F.xxy(iy,ix);
        fxxt = F.xxt(iy,ix);
        fxyy = F.xyy(iy,ix);
        fyyy = F.yyy(iy,ix);
        fyyt = F.yyt(iy,ix);
        fxtt = F.xtt(iy,ix);
        fytt = F.ytt(iy,ix);
        fttt = F.ttt(iy,ix);
        fxyt = F.xyt(iy,ix);
        % Setup space for the 24 linear independent constraint equations,
        % in particular their right-hand-side (vector B) and left-hand-side
        % (matrix A).
        B = zeros(24,1);
        A = zeros(24,8); % u, ux, uy, ut, v, vx, vy, vt
% eqNo. 01:  ( +fx -fxx -fxy -fxt +.5*fxxx +fxxy +fxxt +.5*fxyy +fxyt +.5*fxtt) * u 
%          + ( -fx +fxx +fxy +fxt -.5*fxxx -fxxy -fxxt -.5*fxyy -fxyt -.5*fxtt) * ux 
%          + ( -fx +fxx +fxy +fxt -.5*fxxx -fxxy -fxxt -.5*fxyy -fxyt -.5*fxtt) * uy 
%          + ( -fx +fxx +fxy +fxt -.5*fxxx -fxxy -fxxt -.5*fxyy -fxyt -.5*fxtt) * ut 
%          + ( +fy -fxy -fyy -fyt +.5*fxxy +fxyy +fxyt +.5*fyyy +fyyt +.5*fytt) * v 
%          + ( -fy +fxy +fyy +fyt -.5*fxxy -fxyy -fxyt -.5*fyyy -fyyt -.5*fytt) * vx 
%          + ( -fy +fxy +fyy +fyt -.5*fxxy -fxyy -fxyt -.5*fyyy -fyyt -.5*fytt) * vy 
%          + ( -fy +fxy +fyy +fyt -.5*fxxy -fxyy -fxyt -.5*fyyy -fyyt -.5*fytt) * vt 
%          + ( +ft -fxt -fyt -ftt +.5*fxxt +fxyt +fxtt +.5*fyyt +fytt +.5*fttt) * -1
        A(1,1) = +fx -fxx -fxy -fxt +.5*fxxx +fxxy +fxxt +.5*fxyy +fxyt +.5*fxtt;
        A(1,2) = -fx +fxx +fxy +fxt -.5*fxxx -fxxy -fxxt -.5*fxyy -fxyt -.5*fxtt;
        A(1,3) = -fx +fxx +fxy +fxt -.5*fxxx -fxxy -fxxt -.5*fxyy -fxyt -.5*fxtt;
        A(1,4) = -fx +fxx +fxy +fxt -.5*fxxx -fxxy -fxxt -.5*fxyy -fxyt -.5*fxtt;
        A(1,5) = +fy -fxy -fyy -fyt +.5*fxxy +fxyy +fxyt +.5*fyyy +fyyt +.5*fytt;
        A(1,6) = -fy +fxy +fyy +fyt -.5*fxxy -fxyy -fxyt -.5*fyyy -fyyt -.5*fytt;
        A(1,7) = -fy +fxy +fyy +fyt -.5*fxxy -fxyy -fxyt -.5*fyyy -fyyt -.5*fytt;
        A(1,8) = -fy +fxy +fyy +fyt -.5*fxxy -fxyy -fxyt -.5*fyyy -fyyt -.5*fytt;
        B(1)   = +ft -fxt -fyt -ftt +.5*fxxt +fxyt +fxtt +.5*fyyt +fytt +.5*fttt;
% eqNo. 03:  ( +fx +fxx -fxy -fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt) * u 
%          + ( +fx +fxx -fxy -fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt) * ux 
%          + ( -fx -fxx +fxy +fxt -.5*fxxx +fxxy +fxxt -.5*fxyy -fxyt -.5*fxtt) * uy 
%          + ( -fx -fxx +fxy +fxt -.5*fxxx +fxxy +fxxt -.5*fxyy -fxyt -.5*fxtt) * ut 
%          + ( +fy +fxy -fyy -fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt) * v 
%          + ( +fy +fxy -fyy -fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt) * vx 
%          + ( -fy -fxy +fyy +fyt -.5*fxxy +fxyy +fxyt -.5*fyyy -fyyt -.5*fytt) * vy 
%          + ( -fy -fxy +fyy +fyt -.5*fxxy +fxyy +fxyt -.5*fyyy -fyyt -.5*fytt) * vt 
%          + ( +ft +fxt -fyt -ftt +.5*fxxt -fxyt -fxtt +.5*fyyt +fytt +.5*fttt) * -1
        A(2,1) = +fx +fxx -fxy -fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt;
        A(2,2) = +fx +fxx -fxy -fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt;
        A(2,3) = -fx -fxx +fxy +fxt -.5*fxxx +fxxy +fxxt -.5*fxyy -fxyt -.5*fxtt;
        A(2,4) = -fx -fxx +fxy +fxt -.5*fxxx +fxxy +fxxt -.5*fxyy -fxyt -.5*fxtt;
        A(2,5) = +fy +fxy -fyy -fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt;
        A(2,6) = +fy +fxy -fyy -fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt;
        A(2,7) = -fy -fxy +fyy +fyt -.5*fxxy +fxyy +fxyt -.5*fyyy -fyyt -.5*fytt;
        A(2,8) = -fy -fxy +fyy +fyt -.5*fxxy +fxyy +fxyt -.5*fyyy -fyyt -.5*fytt;
        B(2)   = +ft +fxt -fyt -ftt +.5*fxxt -fxyt -fxtt +.5*fyyt +fytt +.5*fttt;      
% eqNo. 07:  ( +fx -fxx +fxy -fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt) * u 
%          + ( -fx +fxx -fxy +fxt -.5*fxxx +fxxy -fxxt -.5*fxyy +fxyt -.5*fxtt) * ux 
%          + ( +fx -fxx +fxy -fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt) * uy 
%          + ( -fx +fxx -fxy +fxt -.5*fxxx +fxxy -fxxt -.5*fxyy +fxyt -.5*fxtt) * ut 
%          + ( +fy -fxy +fyy -fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt) * v 
%          + ( -fy +fxy -fyy +fyt -.5*fxxy +fxyy -fxyt -.5*fyyy +fyyt -.5*fytt) * vx 
%          + ( +fy -fxy +fyy -fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt) * vy 
%          + ( -fy +fxy -fyy +fyt -.5*fxxy +fxyy -fxyt -.5*fyyy +fyyt -.5*fytt) * vt 
%          + ( +ft -fxt +fyt -ftt +.5*fxxt -fxyt +fxtt +.5*fyyt -fytt +.5*fttt) * -1
        A(3,1) = +fx -fxx +fxy -fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt;
        A(3,2) = -fx +fxx -fxy +fxt -.5*fxxx +fxxy -fxxt -.5*fxyy +fxyt -.5*fxtt;
        A(3,3) = +fx -fxx +fxy -fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt;
        A(3,4) = -fx +fxx -fxy +fxt -.5*fxxx +fxxy -fxxt -.5*fxyy +fxyt -.5*fxtt;
        A(3,5) = +fy -fxy +fyy -fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt;
        A(3,6) = -fy +fxy -fyy +fyt -.5*fxxy +fxyy -fxyt -.5*fyyy +fyyt -.5*fytt;
        A(3,7) = +fy -fxy +fyy -fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt;
        A(3,8) = -fy +fxy -fyy +fyt -.5*fxxy +fxyy -fxyt -.5*fyyy +fyyt -.5*fytt;
        B(3)   = +ft -fxt +fyt -ftt +.5*fxxt -fxyt +fxtt +.5*fyyt -fytt +.5*fttt;
% eqNo. 09:  ( +fx +fxx +fxy -fxt +.5*fxxx +fxxy -fxxt +.5*fxyy -fxyt +.5*fxtt) * u 
%          + ( +fx +fxx +fxy -fxt +.5*fxxx +fxxy -fxxt +.5*fxyy -fxyt +.5*fxtt) * ux 
%          + ( +fx +fxx +fxy -fxt +.5*fxxx +fxxy -fxxt +.5*fxyy -fxyt +.5*fxtt) * uy 
%          + ( -fx -fxx -fxy +fxt -.5*fxxx -fxxy +fxxt -.5*fxyy +fxyt -.5*fxtt) * ut 
%          + ( +fy +fxy +fyy -fyt +.5*fxxy +fxyy -fxyt +.5*fyyy -fyyt +.5*fytt) * v 
%          + ( +fy +fxy +fyy -fyt +.5*fxxy +fxyy -fxyt +.5*fyyy -fyyt +.5*fytt) * vx 
%          + ( +fy +fxy +fyy -fyt +.5*fxxy +fxyy -fxyt +.5*fyyy -fyyt +.5*fytt) * vy 
%          + ( -fy -fxy -fyy +fyt -.5*fxxy -fxyy +fxyt -.5*fyyy +fyyt -.5*fytt) * vt 
%          + ( +ft +fxt +fyt -ftt +.5*fxxt +fxyt -fxtt +.5*fyyt -fytt +.5*fttt) * -1 
% kick out one, e.g. 9.
% eqNo. 14:  ( +fx) * u 
%          + () * ux 
%          + () * uy 
%          + () * ut 
%          + ( +fy) * v 
%          + () * vx 
%          + () * vy 
%          + () * vt 
%          + ( +ft) * -1 + 
        A(4,1) = +fx;
        A(4,5) = +fy;
        B(4)   = +ft;
% eqNo. 19:  ( +fx -fxx -fxy +fxt +.5*fxxx +fxxy -fxxt +.5*fxyy -fxyt +.5*fxtt) * u 
%          + ( -fx +fxx +fxy -fxt -.5*fxxx -fxxy +fxxt -.5*fxyy +fxyt -.5*fxtt) * ux 
%          + ( -fx +fxx +fxy -fxt -.5*fxxx -fxxy +fxxt -.5*fxyy +fxyt -.5*fxtt) * uy 
%          + ( +fx -fxx -fxy +fxt +.5*fxxx +fxxy -fxxt +.5*fxyy -fxyt +.5*fxtt) * ut 
%          + ( +fy -fxy -fyy +fyt +.5*fxxy +fxyy -fxyt +.5*fyyy -fyyt +.5*fytt) * v 
%          + ( -fy +fxy +fyy -fyt -.5*fxxy -fxyy +fxyt -.5*fyyy +fyyt -.5*fytt) * vx 
%          + ( -fy +fxy +fyy -fyt -.5*fxxy -fxyy +fxyt -.5*fyyy +fyyt -.5*fytt) * vy 
%          + ( +fy -fxy -fyy +fyt +.5*fxxy +fxyy -fxyt +.5*fyyy -fyyt +.5*fytt) * vt 
%          + ( +ft -fxt -fyt +ftt +.5*fxxt +fxyt -fxtt +.5*fyyt -fytt +.5*fttt) * -1 + 
        A(5,1) = +fx -fxx -fxy +fxt +.5*fxxx +fxxy -fxxt +.5*fxyy -fxyt +.5*fxtt;
        A(5,2) = -fx +fxx +fxy -fxt -.5*fxxx -fxxy +fxxt -.5*fxyy +fxyt -.5*fxtt;
        A(5,3) = -fx +fxx +fxy -fxt -.5*fxxx -fxxy +fxxt -.5*fxyy +fxyt -.5*fxtt;
        A(5,4) = +fx -fxx -fxy +fxt +.5*fxxx +fxxy -fxxt +.5*fxyy -fxyt +.5*fxtt;
        A(5,5) = +fy -fxy -fyy +fyt +.5*fxxy +fxyy -fxyt +.5*fyyy -fyyt +.5*fytt;
        A(5,6) = -fy +fxy +fyy -fyt -.5*fxxy -fxyy +fxyt -.5*fyyy +fyyt -.5*fytt;
        A(5,7) = -fy +fxy +fyy -fyt -.5*fxxy -fxyy +fxyt -.5*fyyy +fyyt -.5*fytt;
        A(5,8) = +fy -fxy -fyy +fyt +.5*fxxy +fxyy -fxyt +.5*fyyy -fyyt +.5*fytt;
        B(5)   = +ft -fxt -fyt +ftt +.5*fxxt +fxyt -fxtt +.5*fyyt -fytt +.5*fttt;
% eqNo. 21:  ( +fx +fxx -fxy +fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt) * u 
%          + ( +fx +fxx -fxy +fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt) * ux 
%          + ( -fx -fxx +fxy -fxt -.5*fxxx +fxxy -fxxt -.5*fxyy +fxyt -.5*fxtt) * uy 
%          + ( +fx +fxx -fxy +fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt) * ut 
%          + ( +fy +fxy -fyy +fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt) * v 
%          + ( +fy +fxy -fyy +fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt) * vx 
%          + ( -fy -fxy +fyy -fyt -.5*fxxy +fxyy -fxyt -.5*fyyy +fyyt -.5*fytt) * vy 
%          + ( +fy +fxy -fyy +fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt) * vt 
%          + ( +ft +fxt -fyt +ftt +.5*fxxt -fxyt +fxtt +.5*fyyt -fytt +.5*fttt) * -1
        A(6,1) = +fx +fxx -fxy +fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt;
        A(6,2) = +fx +fxx -fxy +fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt;
        A(6,3) = -fx -fxx +fxy -fxt -.5*fxxx +fxxy -fxxt -.5*fxyy +fxyt -.5*fxtt;
        A(6,4) = +fx +fxx -fxy +fxt +.5*fxxx -fxxy +fxxt +.5*fxyy -fxyt +.5*fxtt;
        A(6,5) = +fy +fxy -fyy +fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt;
        A(6,6) = +fy +fxy -fyy +fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt;
        A(6,7) = -fy -fxy +fyy -fyt -.5*fxxy +fxyy -fxyt -.5*fyyy +fyyt -.5*fytt;
        A(6,8) = +fy +fxy -fyy +fyt +.5*fxxy -fxyy +fxyt +.5*fyyy -fyyt +.5*fytt;
        B(6)   = +ft +fxt -fyt +ftt +.5*fxxt -fxyt +fxtt +.5*fyyt -fytt +.5*fttt;
% eqNo. 25:  ( +fx -fxx +fxy +fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt) * u 
%          + ( -fx +fxx -fxy -fxt -.5*fxxx +fxxy +fxxt -.5*fxyy -fxyt -.5*fxtt) * ux 
%          + ( +fx -fxx +fxy +fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt) * uy 
%          + ( +fx -fxx +fxy +fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt) * ut 
%          + ( +fy -fxy +fyy +fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt) * v 
%          + ( -fy +fxy -fyy -fyt -.5*fxxy +fxyy +fxyt -.5*fyyy -fyyt -.5*fytt) * vx 
%          + ( +fy -fxy +fyy +fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt) * vy 
%          + ( +fy -fxy +fyy +fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt) * vt 
%          + ( +ft -fxt +fyt +ftt +.5*fxxt -fxyt -fxtt +.5*fyyt +fytt +.5*fttt) * -1
        A(7,1) = +fx -fxx +fxy +fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt;
        A(7,2) = -fx +fxx -fxy -fxt -.5*fxxx +fxxy +fxxt -.5*fxyy -fxyt -.5*fxtt;
        A(7,3) = +fx -fxx +fxy +fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt;
        A(7,4) = +fx -fxx +fxy +fxt +.5*fxxx -fxxy -fxxt +.5*fxyy +fxyt +.5*fxtt;
        A(7,5) = +fy -fxy +fyy +fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt;
        A(7,6) = -fy +fxy -fyy -fyt -.5*fxxy +fxyy +fxyt -.5*fyyy -fyyt -.5*fytt;
        A(7,7) = +fy -fxy +fyy +fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt;
        A(7,8) = +fy -fxy +fyy +fyt +.5*fxxy -fxyy -fxyt +.5*fyyy +fyyt +.5*fytt;
        B(7)   = +ft -fxt +fyt +ftt +.5*fxxt -fxyt -fxtt +.5*fyyt +fytt +.5*fttt;
% eqNo. 27:  ( +fx +fxx +fxy +fxt +.5*fxxx +fxxy +fxxt +.5*fxyy +fxyt +.5*fxtt) * u 
%          + ( +fx +fxx +fxy +fxt +.5*fxxx +fxxy +fxxt +.5*fxyy +fxyt +.5*fxtt) * ux 
%          + ( +fx +fxx +fxy +fxt +.5*fxxx +fxxy +fxxt +.5*fxyy +fxyt +.5*fxtt) * uy 
%          + ( +fx +fxx +fxy +fxt +.5*fxxx +fxxy +fxxt +.5*fxyy +fxyt +.5*fxtt) * ut 
%          + ( +fy +fxy +fyy +fyt +.5*fxxy +fxyy +fxyt +.5*fyyy +fyyt +.5*fytt) * v 
%          + ( +fy +fxy +fyy +fyt +.5*fxxy +fxyy +fxyt +.5*fyyy +fyyt +.5*fytt) * vx 
%          + ( +fy +fxy +fyy +fyt +.5*fxxy +fxyy +fxyt +.5*fyyy +fyyt +.5*fytt) * vy 
%          + ( +fy +fxy +fyy +fyt +.5*fxxy +fxyy +fxyt +.5*fyyy +fyyt +.5*fytt) * vt 
%          + ( +ft +fxt +fyt +ftt +.5*fxxt +fxyt +fxtt +.5*fyyt +fytt +.5*fttt) * -1
% kick out one, e.g. 27.
% eqNo. 02:  ( +fx -fxy -fxt +.5*fxyy +fxyt +.5*fxtt) * u 
%          + () * ux 
%          + ( -fx +fxy +fxt -.5*fxyy -fxyt -.5*fxtt) * uy 
%          + ( -fx +fxy +fxt -.5*fxyy -fxyt -.5*fxtt) * ut 
%          + ( +fy -fyy -fyt +.5*fyyy +fyyt +.5*fytt) * v 
%          + () * vx 
%          + ( -fy +fyy +fyt -.5*fyyy -fyyt -.5*fytt) * vy 
%          + ( -fy +fyy +fyt -.5*fyyy -fyyt -.5*fytt) * vt 
%          + ( +ft -fyt -ftt +.5*fyyt +fytt +.5*fttt) * -1
        A(8,1) = +fx -fxy -fxt +.5*fxyy +fxyt +.5*fxtt;
        A(8,3) = -fx +fxy +fxt -.5*fxyy -fxyt -.5*fxtt;
        A(8,4) = -fx +fxy +fxt -.5*fxyy -fxyt -.5*fxtt;
        A(8,5) = +fy -fyy -fyt +.5*fyyy +fyyt +.5*fytt;
        A(8,7) = -fy +fyy +fyt -.5*fyyy -fyyt -.5*fytt;
        A(8,8) = -fy +fyy +fyt -.5*fyyy -fyyt -.5*fytt;
        B(8)   = +ft -fyt -ftt +.5*fyyt +fytt +.5*fttt;
% eqNo. 08:  ( +fx +fxy -fxt +.5*fxyy -fxyt +.5*fxtt) * u 
%          + () * ux 
%          + ( +fx +fxy -fxt +.5*fxyy -fxyt +.5*fxtt) * uy 
%          + ( -fx -fxy +fxt -.5*fxyy +fxyt -.5*fxtt) * ut 
%          + ( +fy +fyy -fyt +.5*fyyy -fyyt +.5*fytt) * v 
%          + () * vx 
%          + ( +fy +fyy -fyt +.5*fyyy -fyyt +.5*fytt) * vy 
%          + ( -fy -fyy +fyt -.5*fyyy +fyyt -.5*fytt) * vt 
%          + ( +ft +fyt -ftt +.5*fyyt -fytt +.5*fttt) * -1
        A(9,1) = +fx +fxy -fxt +.5*fxyy -fxyt +.5*fxtt;
        A(9,3) = +fx +fxy -fxt +.5*fxyy -fxyt +.5*fxtt;
        A(9,4) = -fx -fxy +fxt -.5*fxyy +fxyt -.5*fxtt;
        A(9,5) = +fy +fyy -fyt +.5*fyyy -fyyt +.5*fytt;
        A(9,7) = +fy +fyy -fyt +.5*fyyy -fyyt +.5*fytt;
        A(9,8) = -fy -fyy +fyt -.5*fyyy +fyyt -.5*fytt;
        B(9)   = +ft +fyt -ftt +.5*fyyt -fytt +.5*fttt;
% eqNo. 04:  ( +fx -fxx -fxt +.5*fxxx +fxxt +.5*fxtt) * u 
%          + ( -fx +fxx +fxt -.5*fxxx -fxxt -.5*fxtt) * ux 
%          + () * uy 
%          + ( -fx +fxx +fxt -.5*fxxx -fxxt -.5*fxtt) * ut 
%          + ( +fy -fxy -fyt +.5*fxxy +fxyt +.5*fytt) * v 
%          + ( -fy +fxy +fyt -.5*fxxy -fxyt -.5*fytt) * vx 
%          + () * vy 
%          + ( -fy +fxy +fyt -.5*fxxy -fxyt -.5*fytt) * vt 
%          + ( +ft -fxt -ftt +.5*fxxt +fxtt +.5*fttt) * -1
        A(10,1) = +fx -fxx -fxt +.5*fxxx +fxxt +.5*fxtt;
        A(10,2) = -fx +fxx +fxt -.5*fxxx -fxxt -.5*fxtt;
        A(10,4) = -fx +fxx +fxt -.5*fxxx -fxxt -.5*fxtt;
        A(10,5) = +fy -fxy -fyt +.5*fxxy +fxyt +.5*fytt;
        A(10,6) = -fy +fxy +fyt -.5*fxxy -fxyt -.5*fytt;
        A(10,8) = -fy +fxy +fyt -.5*fxxy -fxyt -.5*fytt;
        B(10)   = +ft -fxt -ftt +.5*fxxt +fxtt +.5*fttt;
% eqNo. 06:  ( +fx +fxx -fxt +.5*fxxx -fxxt +.5*fxtt) * u 
%          + ( +fx +fxx -fxt +.5*fxxx -fxxt +.5*fxtt) * ux 
%          + () * uy 
%          + ( -fx -fxx +fxt -.5*fxxx +fxxt -.5*fxtt) * ut 
%          + ( +fy +fxy -fyt +.5*fxxy -fxyt +.5*fytt) * v 
%          + ( +fy +fxy -fyt +.5*fxxy -fxyt +.5*fytt) * vx 
%          + () * vy 
%          + ( -fy -fxy +fyt -.5*fxxy +fxyt -.5*fytt) * vt 
%          + ( +ft +fxt -ftt +.5*fxxt -fxtt +.5*fttt) * -1
        A(11,1) = +fx +fxx -fxt +.5*fxxx -fxxt +.5*fxtt;
        A(11,2) = +fx +fxx -fxt +.5*fxxx -fxxt +.5*fxtt;
        A(11,4) = -fx -fxx +fxt -.5*fxxx +fxxt -.5*fxtt;
        A(11,5) = +fy +fxy -fyt +.5*fxxy -fxyt +.5*fytt;
        A(11,6) = +fy +fxy -fyt +.5*fxxy -fxyt +.5*fytt;
        A(11,8) = -fy -fxy +fyt -.5*fxxy +fxyt -.5*fytt;
        B(11)   = +ft +fxt -ftt +.5*fxxt -fxtt +.5*fttt;
% eqNo. 10:  ( +fx -fxx -fxy +.5*fxxx +fxxy +.5*fxyy) * u 
%          + ( -fx +fxx +fxy -.5*fxxx -fxxy -.5*fxyy) * ux 
%          + ( -fx +fxx +fxy -.5*fxxx -fxxy -.5*fxyy) * uy 
%          + () * ut 
%          + ( +fy -fxy -fyy +.5*fxxy +fxyy +.5*fyyy) * v 
%          + ( -fy +fxy +fyy -.5*fxxy -fxyy -.5*fyyy) * vx 
%          + ( -fy +fxy +fyy -.5*fxxy -fxyy -.5*fyyy) * vy 
%          + () * vt 
%          + ( +ft -fxt -fyt +.5*fxxt +fxyt +.5*fyyt) * -1
        A(12,1) = +fx -fxx -fxy +.5*fxxx +fxxy +.5*fxyy;
        A(12,2) = +fx +fxx -fxt +.5*fxxx -fxxt +.5*fxtt;
        A(12,3) = -fx +fxx +fxy -.5*fxxx -fxxy -.5*fxyy;
        A(12,5) = +fy -fxy -fyy +.5*fxxy +fxyy +.5*fyyy;
        A(12,6) = -fy +fxy +fyy -.5*fxxy -fxyy -.5*fyyy;
        A(12,7) = -fy +fxy +fyy -.5*fxxy -fxyy -.5*fyyy;
        B(12)   = +ft -fxt -fyt +.5*fxxt +fxyt +.5*fyyt;
% eqNo. 12:  ( +fx +fxx -fxy +.5*fxxx -fxxy +.5*fxyy) * u 
%          + ( +fx +fxx -fxy +.5*fxxx -fxxy +.5*fxyy) * ux 
%          + ( -fx -fxx +fxy -.5*fxxx +fxxy -.5*fxyy) * uy 
%          + () * ut 
%          + ( +fy +fxy -fyy +.5*fxxy -fxyy +.5*fyyy) * v 
%          + ( +fy +fxy -fyy +.5*fxxy -fxyy +.5*fyyy) * vx 
%          + ( -fy -fxy +fyy -.5*fxxy +fxyy -.5*fyyy) * vy 
%          + () * vt 
%          + ( +ft +fxt -fyt +.5*fxxt -fxyt +.5*fyyt) * -1
        A(13,1) = +fx +fxx -fxy +.5*fxxx -fxxy +.5*fxyy;
        A(13,2) = +fx +fxx -fxy +.5*fxxx -fxxy +.5*fxyy;
        A(13,3) = -fx -fxx +fxy -.5*fxxx +fxxy -.5*fxyy;
        A(13,5) = +fy +fxy -fyy +.5*fxxy -fxyy +.5*fyyy;
        A(13,6) = +fy +fxy -fyy +.5*fxxy -fxyy +.5*fyyy;
        A(13,7) = -fy -fxy +fyy -.5*fxxy +fxyy -.5*fyyy;
        B(13)   = +ft +fxt -fyt +.5*fxxt -fxyt +.5*fyyt;
% eqNo. 16:  ( +fx -fxx +fxy +.5*fxxx -fxxy +.5*fxyy) * u 
%          + ( -fx +fxx -fxy -.5*fxxx +fxxy -.5*fxyy) * ux 
%          + ( +fx -fxx +fxy +.5*fxxx -fxxy +.5*fxyy) * uy 
%          + () * ut 
%          + ( +fy -fxy +fyy +.5*fxxy -fxyy +.5*fyyy) * v 
%          + ( -fy +fxy -fyy -.5*fxxy +fxyy -.5*fyyy) * vx 
%          + ( +fy -fxy +fyy +.5*fxxy -fxyy +.5*fyyy) * vy 
%          + () * vt 
%          + ( +ft -fxt +fyt +.5*fxxt -fxyt +.5*fyyt) * -1
        A(14,1) = +fx -fxx +fxy +.5*fxxx -fxxy +.5*fxyy;
        A(14,2) = -fx +fxx -fxy -.5*fxxx +fxxy -.5*fxyy;
        A(14,3) = +fx -fxx +fxy +.5*fxxx -fxxy +.5*fxyy;
        A(14,5) = +fy -fxy +fyy +.5*fxxy -fxyy +.5*fyyy;
        A(14,6) = -fy +fxy -fyy -.5*fxxy +fxyy -.5*fyyy;
        A(14,7) = +fy -fxy +fyy +.5*fxxy -fxyy +.5*fyyy;
        B(14)   = +ft -fxt +fyt +.5*fxxt -fxyt +.5*fyyt;
% eqNo. 18:  ( +fx +fxx +fxy +.5*fxxx +fxxy +.5*fxyy) * u 
%          + ( +fx +fxx +fxy +.5*fxxx +fxxy +.5*fxyy) * ux 
%          + ( +fx +fxx +fxy +.5*fxxx +fxxy +.5*fxyy) * uy 
%          + () * ut 
%          + ( +fy +fxy +fyy +.5*fxxy +fxyy +.5*fyyy) * v 
%          + ( +fy +fxy +fyy +.5*fxxy +fxyy +.5*fyyy) * vx 
%          + ( +fy +fxy +fyy +.5*fxxy +fxyy +.5*fyyy) * vy 
%          + () * vt 
%          + ( +ft +fxt +fyt +.5*fxxt +fxyt +.5*fyyt) * -1
% kick out one, e.g. 18.
% eqNo. 22:  ( +fx -fxx +fxt +.5*fxxx -fxxt +.5*fxtt) * u 
%          + ( -fx +fxx -fxt -.5*fxxx +fxxt -.5*fxtt) * ux 
%          + () * uy 
%          + ( +fx -fxx +fxt +.5*fxxx -fxxt +.5*fxtt) * ut 
%          + ( +fy -fxy +fyt +.5*fxxy -fxyt +.5*fytt) * v 
%          + ( -fy +fxy -fyt -.5*fxxy +fxyt -.5*fytt) * vx 
%          + () * vy 
%          + ( +fy -fxy +fyt +.5*fxxy -fxyt +.5*fytt) * vt 
%          + ( +ft -fxt +ftt +.5*fxxt -fxtt +.5*fttt) * -1
        A(15,1) = +fx -fxx +fxt +.5*fxxx -fxxt +.5*fxtt;
        A(15,2) = -fx +fxx -fxt -.5*fxxx +fxxt -.5*fxtt;
        A(15,4) = +fx -fxx +fxt +.5*fxxx -fxxt +.5*fxtt;
        A(15,5) = +fy -fxy +fyt +.5*fxxy -fxyt +.5*fytt;
        A(15,6) = -fy +fxy -fyt -.5*fxxy +fxyt -.5*fytt;
        A(15,8) = +fy -fxy +fyt +.5*fxxy -fxyt +.5*fytt;
        B(15)   = +ft -fxt +ftt +.5*fxxt -fxtt +.5*fttt;
% eqNo. 24:  ( +fx +fxx +fxt +.5*fxxx +fxxt +.5*fxtt) * u 
%          + ( +fx +fxx +fxt +.5*fxxx +fxxt +.5*fxtt) * ux 
%          + () * uy 
%          + ( +fx +fxx +fxt +.5*fxxx +fxxt +.5*fxtt) * ut 
%          + ( +fy +fxy +fyt +.5*fxxy +fxyt +.5*fytt) * v 
%          + ( +fy +fxy +fyt +.5*fxxy +fxyt +.5*fytt) * vx 
%          + () * vy 
%          + ( +fy +fxy +fyt +.5*fxxy +fxyt +.5*fytt) * vt 
%          + ( +ft +fxt +ftt +.5*fxxt +fxtt +.5*fttt) * -1
        A(16,1) = +fx +fxx +fxt +.5*fxxx +fxxt +.5*fxtt;
        A(16,2) = +fx +fxx +fxt +.5*fxxx +fxxt +.5*fxtt;
        A(16,4) = +fx +fxx +fxt +.5*fxxx +fxxt +.5*fxtt;
        A(16,5) = +fy +fxy +fyt +.5*fxxy +fxyt +.5*fytt;
        A(16,6) = +fy +fxy +fyt +.5*fxxy +fxyt +.5*fytt;
        A(16,8) = +fy +fxy +fyt +.5*fxxy +fxyt +.5*fytt;
        B(16)   = +ft +fxt +ftt +.5*fxxt +fxtt +.5*fttt;
% eqNo. 20:  ( +fx -fxy +fxt +.5*fxyy -fxyt +.5*fxtt) * u 
%          + () * ux 
%          + ( -fx +fxy -fxt -.5*fxyy +fxyt -.5*fxtt) * uy 
%          + ( +fx -fxy +fxt +.5*fxyy -fxyt +.5*fxtt) * ut 
%          + ( +fy -fyy +fyt +.5*fyyy -fyyt +.5*fytt) * v 
%          + () * vx 
%          + ( -fy +fyy -fyt -.5*fyyy +fyyt -.5*fytt) * vy 
%          + ( +fy -fyy +fyt +.5*fyyy -fyyt +.5*fytt) * vt 
%          + ( +ft -fyt +ftt +.5*fyyt -fytt +.5*fttt) * -1
        A(17,1) = +fx -fxy +fxt +.5*fxyy -fxyt +.5*fxtt;
        A(17,3) = -fx +fxy -fxt -.5*fxyy +fxyt -.5*fxtt;
        A(17,4) = +fx -fxy +fxt +.5*fxyy -fxyt +.5*fxtt;
        A(17,5) = +fy -fyy +fyt +.5*fyyy -fyyt +.5*fytt;
        A(17,7) = -fy +fyy -fyt -.5*fyyy +fyyt -.5*fytt;
        A(17,8) = +fy -fyy +fyt +.5*fyyy -fyyt +.5*fytt;
        B(17)   = +ft -fyt +ftt +.5*fyyt -fytt +.5*fttt;
% eqNo. 26:  ( +fx +fxy +fxt +.5*fxyy +fxyt +.5*fxtt) * u 
%          + () * ux 
%          + ( +fx +fxy +fxt +.5*fxyy +fxyt +.5*fxtt) * uy 
%          + ( +fx +fxy +fxt +.5*fxyy +fxyt +.5*fxtt) * ut 
%          + ( +fy +fyy +fyt +.5*fyyy +fyyt +.5*fytt) * v 
%          + () * vx 
%          + ( +fy +fyy +fyt +.5*fyyy +fyyt +.5*fytt) * vy 
%          + ( +fy +fyy +fyt +.5*fyyy +fyyt +.5*fytt) * vt 
%          + ( +ft +fyt +ftt +.5*fyyt +fytt +.5*fttt) * -1
        A(18,1) = +fx +fxy +fxt +.5*fxyy +fxyt +.5*fxtt;
        A(18,3) = +fx +fxy +fxt +.5*fxyy +fxyt +.5*fxtt;
        A(18,4) = +fy +fyy +fyt +.5*fyyy +fyyt +.5*fytt;
        A(18,5) = +fy +fyy +fyt +.5*fyyy +fyyt +.5*fytt;
        A(18,7) = +fy +fyy +fyt +.5*fyyy +fyyt +.5*fytt;
        A(18,8) = +fy +fyy +fyt +.5*fyyy +fyyt +.5*fytt;
        B(18)   = +ft +fyt +ftt +.5*fyyt +fytt +.5*fttt;
% eqNo. 05:  ( +fx -fxt +.5*fxtt) * u 
%          + () * ux 
%          + () * uy 
%          + ( -fx +fxt -.5*fxtt) * ut 
%          + ( +fy -fyt +.5*fytt) * v 
%          + () * vx 
%          + () * vy 
%          + ( -fy +fyt -.5*fytt) * vt 
%          + ( +ft -ftt +.5*fttt) * -1
        A(19,1) = +fx -fxt +.5*fxtt;
        A(19,4) = -fx +fxt -.5*fxtt;
        A(19,5) = +fy -fyt +.5*fytt;
        A(19,8) = -fy +fyt -.5*fytt;
        B(19)   = +ft -ftt +.5*fttt;
% eqNo. 23:  ( +fx +fxt +.5*fxtt) * u 
%          + () * ux 
%          + () * uy 
%          + ( +fx +fxt +.5*fxtt) * ut 
%          + ( +fy +fyt +.5*fytt) * v 
%          + () * vx 
%          + () * vy 
%          + ( +fy +fyt +.5*fytt) * vt 
%          + ( +ft +ftt +.5*fttt) * -1
        A(20,1) = +fx +fxt +.5*fxtt;
        A(20,4) = +fx +fxt +.5*fxtt;
        A(20,5) = +fy +fyt +.5*fytt;
        A(20,8) = +fy +fyt +.5*fytt;
        B(20)   = +ft +ftt +.5*fttt;
% eqNo. 13:  ( +fx -fxx +.5*fxxx) * u 
%          + ( -fx +fxx -.5*fxxx) * ux 
%          + () * uy 
%          + () * ut 
%          + ( +fy -fxy +.5*fxxy) * v 
%          + ( -fy +fxy -.5*fxxy) * vx 
%          + () * vy 
%          + () * vt 
%          + ( +ft -fxt +.5*fxxt) * -1
        A(21,1) = +fx -fxx +.5*fxxx;
        A(21,2) = -fx +fxx -.5*fxxx;
        A(21,5) = +fy -fxy +.5*fxxy;
        A(21,6) = -fy +fxy -.5*fxxy;
        B(21)   = +ft -fxt +.5*fxxt;
% eqNo. 15:  ( +fx +fxx +.5*fxxx) * u 
%          + ( +fx +fxx +.5*fxxx) * ux 
%          + () * uy 
%          + () * ut 
%          + ( +fy +fxy +.5*fxxy) * v 
%          + ( +fy +fxy +.5*fxxy) * vx 
%          + () * vy 
%          + () * vt 
%          + ( +ft +fxt +.5*fxxt) * -1
        A(22,1) = +fx +fxx +.5*fxxx;
        A(22,2) = +fx +fxx +.5*fxxx;
        A(22,5) = +fy +fxy +.5*fxxy;
        A(22,6) = +fy +fxy +.5*fxxy;
        B(22)   = +ft +fxt +.5*fxxt;
% eqNo. 11:  ( +fx -fxy +.5*fxyy) * u 
%          + () * ux 
%          + ( -fx +fxy -.5*fxyy) * uy 
%          + () * ut 
%          + ( +fy -fyy +.5*fyyy) * v 
%          + () * vx 
%          + ( -fy +fyy -.5*fyyy) * vy 
%          + () * vt 
%          + ( +ft -fyt +.5*fyyt) * -1
        A(23,1) = +fx -fxy +.5*fxyy;
        A(23,3) = -fx +fxy -.5*fxyy;
        A(23,5) = +fy -fyy +.5*fyyy;
        A(23,7) = -fy +fyy -.5*fyyy;
        B(23)   = +ft -fyt +.5*fyyt;
% eqNo. 17:  ( +fx +fxy +.5*fxyy) * u 
%          + () * ux 
%          + ( +fx +fxy +.5*fxyy) * uy 
%          + () * ut 
%          + ( +fy +fyy +.5*fyyy) * v 
%          + () * vx 
%          + ( +fy +fyy +.5*fyyy) * vy 
%          + () * vt 
%          + ( +ft +fyt +.5*fyyt) * -1
        A(24,1) = +fx +fxy +.5*fxyy;
        A(24,3) = +fx +fxy +.5*fxyy;
        A(24,5) = +fy +fyy +.5*fyyy;
        A(24,7) = +fy +fyy +.5*fyyy;
        B(24)   = +ft +fyt +.5*fyyt;
        % Compute the least squares solution for the flow and its
        % derivatives.
        X = - pinv(A) * B;
        Dx(iy,ix,:) = X(1:4);
        Dy(iy,ix,:) = X(5:8);
    end
end
