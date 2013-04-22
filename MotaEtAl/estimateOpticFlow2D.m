function [Dx Dy] = estimateOpticFlow2D(ImgSeq, opt)
% estimateOpticFlow2D
%   ImgSeq  - Image sequence as a cube with dimensions: 
%             height x width x frames.
%   opt     - Struture with options:
%             * W       - Weight matrix for sums constraints in the local 
%                         spatio-temporal neighborhood to build the 
%                         structure tensor.
%             * thOne   - Threshold for ONE motion component being present.
%             * thTwo   - Threshold for TWO motions being present.
%             * DiffOp  - Differential operator for spatial and temporal
%                         domain.
%
% RETURN
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             Both flow components of the dimensions: height x width x 2.
%             This method computes the flow for the center frame.
%
% DESCRIPTION
%   An implementation of 
%   Mota, C., Stuke, I. and Barth, E. (2001). Analytical solutions for 
%       multiple motions. In Proceedings of IEEE International Conference 
%       on Image Processing (ICIP'01), Thessaloniki, Greece.
%   In this implementation, I estiamte at maximum two motion components.
%   More details and some parameter values are described in
%   Mota, C., Stuke, I., Aach, T., and Barth, E. (2004). Divide-and-conquer
%       strategies for estimating multiple transparent motions. In Jähne,
%       B. et al. (eds.) Complex Motion, LNCS 3417, 66-77.
%
%   Copyright (C) 2013  Florian Raudies, 01/08/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for parameters of the method.
if nargin<2,                opt         = struct();     end % UNITS
if ~isfield(opt,'thOne'),   opt.thOne   = 0.7;          end % -
if ~isfield(opt,'thTwo'),   opt.thTwo   = 0.7;          end % -
if ~isfield(opt,'DiffOp'),  opt.DiffOp  = [-1 0 1]/2;   end % -
if ~isfield(opt,'W'),       opt.W       = fspecial('gaussian', [5 1], 1); end
% Retrieve parameter values.
thOne   = opt.thOne;
thTwo   = opt.thTwo;
DiffOp  = opt.DiffOp;
W       = opt.W;
% Use Sobel operators to calculate the partial derivatives, see 
% Mota et al. (2004).
yNum        = size(ImgSeq,1);
xNum        = size(ImgSeq,2);
tNum        = size(ImgSeq,3);
DiffOp      = reshape(DiffOp,[1 1 3]);
ktNum       = size(DiffOp,3);
Binomi2x2   = convn([1 2 1],[1 2 1]');
DiffT       = convn(Binomi2x2, DiffOp);
DiffX       = shiftdim(DiffT,1);
DiffY       = shiftdim(DiffT,2);
Wx          = reshape(W,[1 5 1]);
Wy          = reshape(W,[5 1 1]);
Wt          = reshape(W,[1 1 5]);
wtNum       = size(Wt,3);
dtNum       = size(DiffT,3);
% Check if the provided sequence contains enough frames.
if tNum<(wtNum+2*dtNum-2), 
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], wtNum+2*dtNum-2, tNum);
end
% Select the center number of wtNum-2*dtNum-2 frames.
center  = ceil(tNum/2);
SelT    = center + ( -(wtNum-1)/2-(dtNum-1) : +(wtNum-1)/2+(dtNum-1) );
ImgSeq  = ImgSeq(:,:,SelT);
tNum    = size(ImgSeq,3);
% Allocate memory for motions.
Dx1 = NaN(yNum,xNum);
Dy1 = NaN(yNum,xNum);
Dx2 = NaN(yNum,xNum);
Dy2 = NaN(yNum,xNum);
% *************************************************************************
% *************************************************************************
% Assume ONE motion component.
% *************************************************************************
% *************************************************************************
mNum = 1;                   % Motion components
dNum = (mNum+1)*(mNum+2)/2; % Number of derivatives
eNum = (dNum^2+dNum)/2;     % Number of elements in the tensor
% I. Calculate first-order derivatives.
D1      = zeros(yNum, xNum, tNum, dNum);
D1(:,:,:,1) = imfilter(ImgSeq, DiffX, 'same', 'replicate');
D1(:,:,:,2) = imfilter(ImgSeq, DiffY, 'same', 'replicate');
D1(:,:,:,3) = imfilter(ImgSeq, DiffT, 'same', 'replicate');
ValidT  = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
D1      = D1(:,:,ValidT,:);
tNum    = size(D1,3);
SelT    = ceil(tNum/2) + (-floor(wtNum/2) : +floor(wtNum/2));
% II. Calculate components of the generalized structure tensor.
T = zeros(yNum, xNum, eNum); % xx, yy, tt, xy, yt, xt
Index = [1 1; 2 2; 3 3; 1 2; 2 3; 1 3]; % Index of tensor components
%      /t1 t4 t6\
% T = | t4 t2 t5 |
%      \t6 t5 t3/
for ie = 1:eNum,
    k = Index(ie,1);
    l = Index(ie,2);
    T(:,:,ie) = imfilter( ...
                  imfilter( ...
                   convn( D1(:,:,SelT,k).*D1(:,:,SelT,l), Wt, 'valid'), ...
                        Wy, 'same', 'replicate'), ...
                            Wx, 'same', 'replicate');
end
% III. Calculate minors.
M = zeros(yNum, xNum, 3, 3);
% First row.
M(:,:,1,1) = T(:,:,2).*T(:,:,3) - T(:,:,5).^2;          % fx*fx
M(:,:,1,2) = T(:,:,4).*T(:,:,3) - T(:,:,5).*T(:,:,6);   % fx*fy
M(:,:,1,3) = T(:,:,4).*T(:,:,5) - T(:,:,2).*T(:,:,6);   % fx*ft
% Second row.
M(:,:,2,1) = M(:,:,1,2);                                % fx*fy
M(:,:,2,2) = T(:,:,1).*T(:,:,3) - T(:,:,6).^2;          % fy*fy
M(:,:,2,3) = T(:,:,1).*T(:,:,5) - T(:,:,4).*T(:,:,6);   % fy*ft
% Third row.
M(:,:,3,1) = M(:,:,1,3);                                % fx*ft
M(:,:,3,2) = M(:,:,2,3);                                % fy*ft
M(:,:,3,3) = T(:,:,1).*T(:,:,2) - T(:,:,4).^2;          % ft*ft
% Calculate determinant using tensor components and minors.
K = squeeze(T(:,:,1)).*squeeze(M(:,:,1,1)) ...
  - squeeze(T(:,:,4)).*squeeze(M(:,:,1,2)) ...
  + squeeze(T(:,:,6)).*squeeze(M(:,:,1,3));
% Mean calculated for diagonal entries of matrix with minors.
S = 1/3 * squeeze(M(:,:,1,1) + M(:,:,2,2) + M(:,:,3,3));
K = K.^(1/dNum);
S = S.^(1/(dNum-1));
% IV. Calculate one motion from minors.
Dx = squeeze(mean(+M(:,:,:,1)./(M(:,:,:,3)+eps),3)); % 3=ft*ft
Dy = squeeze(mean(-M(:,:,:,2)./(M(:,:,:,3)+eps),3));
% V. Check for one motion. 
OneMotion = K < thOne * S;
% Store result in output matrices.
Dx1(OneMotion) = Dx(OneMotion);
Dy1(OneMotion) = Dy(OneMotion);
% *************************************************************************
% *************************************************************************
% Assume TWO motion components.
% *************************************************************************
% *************************************************************************
tNum = size(D1,3);
mNum = 2;                   % Motion components.
dNum = (mNum+1)*(mNum+2)/2; % Number of derivatives.
eNum = (dNum^2+dNum)/2;     % Number of elements in the tensor.
% I. Calculate 2nd order derivatives.
D2 = zeros(yNum, xNum, tNum, dNum);
D2(:,:,:,1) = imfilter(squeeze(D1(:,:,:,1)), DiffX, 'same', 'replicate'); % xx
D2(:,:,:,2) = imfilter(squeeze(D1(:,:,:,2)), DiffY, 'same', 'replicate'); % yy
D2(:,:,:,3) = imfilter(squeeze(D1(:,:,:,3)), DiffT, 'same', 'replicate'); % tt
D2(:,:,:,4) = imfilter(squeeze(D1(:,:,:,1)), DiffY, 'same', 'replicate'); % xy
D2(:,:,:,5) = imfilter(squeeze(D1(:,:,:,2)), DiffT, 'same', 'replicate'); % yt
D2(:,:,:,6) = imfilter(squeeze(D1(:,:,:,1)), DiffT, 'same', 'replicate'); % xt
ValidT  = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
D2 = D2(:,:,ValidT,:);
% II. Calculate components of the generalized structure tensor.
T = zeros(yNum, xNum, eNum);
%          1    2    3    4    5    6    7    8    9   10   11
Index = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 1 2; 2 3; 3 4; 4 5; 5 6; ...
         1 3; 2 4; 3 5; 4 6; 1 4; 2 5; 3 6; 1 5; 2 6; 1 6];
tNum  = size(D2,3);     
SelT  = ceil(tNum/2) + (-floor(wtNum/2) : +floor(wtNum/2));
for ie = 1:eNum,
    k = Index(ie,1);
    l = Index(ie,2);
    T(:,:,ie) = imfilter( ...
                  imfilter( ...
                   convn( D2(:,:,SelT,k).*D2(:,:,SelT,l), Wt, 'valid'), ...
                        Wy, 'same', 'replicate'), ...
                            Wx, 'same', 'replicate');
end
% 1=fxx*fxx, 2=fyy*fyy, 3=ftt*ftt, 4=fxy*fxy, 5=fyt*fyt, 6=fxt*fxt, 
% 7=fxx*fyy, 8=fyy*ftt, 9=ftt*fxy, 10=fxy*fyt, 11=fyt*fxt, 12=fxx*ftt,
% 13=fyy*fxy, 14=ftt*fyt, 15=fxy*fxt, 16=fxx*fxy, 17=fyy*fyt, 18=ftt*fxt,
% 19=fxx*fyt, 20=fyy*fxt, 21=fxx*fxt.
Index = [ 1  7 12 16 19 21; ...
          7  2  8 13 17 20; ...
         12  8  3  9 14 18; ...
         16 13  9  4 10 15; ...
         19 17 14 10  5 11; ...
         21 20 18 15 11  6];
% III. Calculate minors.
K       = zeros(yNum, xNum);
S       = zeros(yNum, xNum);
Vsum    = zeros(yNum, xNum, dNum); % xx, yy, tt, xy, yt, xt
% For each location.
for iy = 1:yNum,
    for ix = 1:xNum,        
        % Expand the tensor to all components.
        Txy         = squeeze(T(iy,ix,:));
        J           = Txy(Index);
        K(iy,ix)    = det(J);
        M           = minors(J);
        S(iy,ix)    = mean(diag(M));
        % Calculate weights from minors.
        Alpha = M(3,:)/(sum(M(3,:).^2)+eps); % 3=ftt*ftt
        % Calculate Vsum expressions.
        Vsum(iy,ix,:) = sum(repmat(Alpha,[dNum 1])' ...
            .*(repmat( (-1).^(0:dNum-1), [dNum 1]).*M),1);
    end
end
K = K.^(1/dNum);
S = S.^(1/(dNum-1));
% IV. Calculate two motions by solving a complex quadratic polynomial.
% A0 = cxx - cyy + j*cxy
A0 = +squeeze( Vsum(:,:,1) - Vsum(:,:,2) + 1j*Vsum(:,:,4) );
% A1 = -(cxt + jcyt), note that I added '-' to the definition of A1.
A1 = -squeeze( Vsum(:,:,6) + 1j*Vsum(:,:,5) );
% Z_1,2 = -(A1)/2 +- sqrt( ((A1)/2).^2 - A0 )
Disc = sqrt( (A1/2).^2 - A0 );
Z1 = -A1/2 + Disc; 
Z2 = -A1/2 - Disc;
% V. Check for two motions.
TwoMotions = (K < thTwo * S) & ~OneMotion;
% Store motions in output the structure.
Dx1(TwoMotions) = real(Z1(TwoMotions));
Dy1(TwoMotions) = imag(Z1(TwoMotions));
Dx2(TwoMotions) = real(Z2(TwoMotions));
Dy2(TwoMotions) = imag(Z2(TwoMotions));
% Merge matrices by adding a third dimension.
Dx = cat(3, Dx1, Dx2);
Dy = cat(3, Dy1, Dy2);
