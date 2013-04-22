function [Dx Dy] = estimateOpticFlow2D(ImgSeq, opt)
% estimateOpticFlow2D
%   ImgSeq  - Image sequence as a cube with dimensions: 
%             height x width x frames.
%   opt     - Struture with options:
%             * W        - Weight matrix for sums constraints in the local 
%                          neighborhood to build the structure tensor.
%             * DiffX    - Kernel that approximtes the computation of the
%                          partial derivative in x.
%             * DiffY    - Ditto for y.
%             * DiffT    - Ditto for t.
%             * eta      - Regulates the smoothness constraint. A larger 
%                          eta results in an increased smoothness.
%             * omega    - Parameter for the successive over relaxation 
%                          (SOR) method, which is used to solve the sparse 
%                          linear equation system. 0 <= omega <= 2.
%             * iNum     - Number of iterations for the SOR method.
%             * showFlow - If showFlow=1, then the flow per iteration is
%                          displayed.
%
% RETURN
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             Computes the optic flow for all frames of the sequence,
%             thus, Dx and Dy have dimensions: height x width x frames-1.
%
% DESCRIPTION
%   A modern implementation of the idea of 
%   Horn, B.K.P. and Schunck, B.G. (1981). Determining optical flow. 
%       Artificial Intelligence 17, 185-203.
%   The discretization and implementation follows:
%   Bruhn, A., Weickert, J. Kohlberger, T., and Schnörr, C. (2006). A 
%       multigrid platform for real-time motion computation with 
%       discontinuity preserving variational methods. International 
%       Journal of Computer Vision 70(3), 257-277.
%
%   Copyright (C) 2013  Florian Raudies, 01/02/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for paraemters of the method.
if nargin<2,                opt         = struct();                     end
if ~isfield(opt,'W'),       opt.W       = fspecial('gaussian',[9 9],2); end
if ~isfield(opt,'DiffX'),   opt.DiffX   = [-1 8 0 -8 1]/12;             end
if ~isfield(opt,'DiffY'),   opt.DiffY   = opt.DiffX';                   end
if ~isfield(opt,'DiffT'),   opt.DiffT   = reshape([-1 1],[1 1 2]);      end
if ~isfield(opt,'eta'),     opt.eta     = 1/255;                        end
if ~isfield(opt,'omega'),   opt.omega   = 1.5;                          end
if ~isfield(opt,'iNum'),    opt.iNum    = 50;                           end
if ~isfield(opt,'showFlow'),opt.showFlow= 1;                            end
% Retreive default parmeters to be saved in their own variables.
W           = opt.W;
DiffX       = opt.DiffX;
DiffY       = opt.DiffY;
DiffT       = opt.DiffT;
eta         = opt.eta;
omega       = opt.omega;
iNum        = opt.iNum;
showFlow    = opt.showFlow;
% Check if the provided sequence contains at least two frames.
[yNum xNum tNum] = size(ImgSeq);
if tNum<2, 
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], 2, tNum);
end
% Compute the partial derivatives in x, y, and t.
ImgSeqDx = imfilter(ImgSeq, DiffX, 'same', 'replicate');
ImgSeqDy = imfilter(ImgSeq, DiffY, 'same', 'replicate');
ImgSeqDt = convn(ImgSeq, DiffT, 'valid');
% Select the spatial and temporal valid part of the partial derivatives.
ktNum    = size(DiffT,3);
ValidT   = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
ImgSeqDx = squeeze(ImgSeqDx(:,:,ValidT));
ImgSeqDy = squeeze(ImgSeqDy(:,:,ValidT));
tNum     = length(ValidT);
Dx       = squeeze(zeros(yNum, xNum, tNum));
Dy       = squeeze(zeros(yNum, xNum, tNum));
% Compute the coefficients A, B, C, D, E, and F of the structure tensor
%     / A B C \
% J = | B D E |
%     \ C E F /
A = imfilter(ImgSeqDx.*ImgSeqDx,W, 'same', 'replicate');
B = imfilter(ImgSeqDx.*ImgSeqDy,W, 'same', 'replicate');
C = imfilter(ImgSeqDx.*ImgSeqDt,W, 'same', 'replicate');
D = imfilter(ImgSeqDy.*ImgSeqDy,W, 'same', 'replicate');
E = imfilter(ImgSeqDy.*ImgSeqDt,W, 'same', 'replicate');
% F = imfilter(ImgSeqDt.*ImgSeqDt,W, 'same', 'replicate'); % not used.
% Repeat boundary values by 1.
A   = repeatBoundary(A, 1);
B   = repeatBoundary(B, 1);
C   = repeatBoundary(C, 1);
D   = repeatBoundary(D, 1);
E   = repeatBoundary(E, 1);
Dx  = repeatBoundary(Dx, 1);
Dy  = repeatBoundary(Dy, 1);
yNum    = yNum + 2;
xNum    = xNum + 2;
IndexY  = 2 : (yNum-1);
IndexX  = 2 : (xNum-1);
% *************************************************************************
% Compute 2D optic flow.
% *************************************************************************
if tNum==1,
    for iter = 1:iNum,
        for iy = IndexY,
            for ix = IndexX,
                Dx(iy,ix) = (1-omega) * Dx(iy,ix) ...
                    + omega * (Dx(iy-1,ix) + Dx(iy+1,ix)...
                              +Dx(iy,ix-1) + Dx(iy,ix+1)...
                              -1/eta*(B(iy,ix)*Dy(iy,ix) + C(iy,ix))) ...
                             /(1/eta*A(iy,ix) + 4);
                Dy(iy,ix) = (1-omega) * Dy(iy,ix) ...
                    + omega * (Dy(iy-1,ix) + Dy(iy+1,ix) ...
                              +Dy(iy,ix-1) + Dy(iy,ix+1) ...
                              -1/eta*(B(iy,ix)*Dx(iy,ix) + E(iy,ix))) ...
                             /(1/eta*D(iy,ix) + 4);
            end
        end
        % Dx = copyBoundary(Dx, kNum);
        % Dy = copyBoundary(Dy, kNum);
        % Use the faster, direct way to copy boundary values.
        Dx(1,:) = Dx(2,:);    Dx(yNum,:) = Dx(yNum-1,:);
        Dx(:,1) = Dx(:,2);    Dx(:,xNum) = Dx(:,xNum-1);
        Dy(1,:) = Dy(2,:);    Dy(yNum,:) = Dy(yNum-1,:);
        Dy(:,1) = Dy(:,2);    Dy(:,xNum) = Dy(:,xNum-1);
        
        % Optionally, display the estimated flow of the current iteration.
        if showFlow, plotFlow(Dx,Dy,iter,iNum); end
    end
% *************************************************************************
% Compute 3D optic flow.
% *************************************************************************
else
    tNum    = tNum + 2;
    IndexT  = 2 : tNum-1;
    for iter = 1:iNum,
        for iy = IndexY,
            for ix = IndexX,
                for it = IndexT,
                    Dx(iy,ix,it) = (1-omega) * Dx(iy,ix,it) ...
                        + omega * (Dx(iy-1,ix,it) + Dx(iy+1,ix,it) ...
                                  +Dx(iy,ix-1,it) + Dx(iy,ix+1,it) ...
                                  +Dx(iy,ix,it-1) + Dx(iy,ix,it+1) ...
                                  -1/eta*(B(iy,ix,it)*Dy(iy,ix,it) + C(iy,ix,it))) ...
                                 /(1/eta*A(iy,ix,it) + 6);
                    Dy(iy,ix,it) = (1-omega) * Dy(iy,ix,it) ...
                        + omega * (Dy(iy-1,ix,it) + Dy(iy+1,ix,it) ...
                                  +Dy(iy,ix-1,it) + Dy(iy,ix+1,it) ...
                                  +Dy(iy,ix,it-1) + Dy(iy,ix,it+1) ...
                                  -1/eta*(B(iy,ix,it)*Dx(iy,ix,it) + E(iy,ix,it))) ...
                                 /(1/eta*D(iy,ix,it) + 6);
                end
            end
        end
        Dx = copyBoundary(Dx, 1);
        Dy = copyBoundary(Dy, 1);
        
        % Optionally, display the estimated flow (1st frame) of the current iteration.
        if showFlow, plotFlow(Dx(IndexY,IndexX,1),Dy(IndexY,IndexX,1),iter,iNum); end    
    end
end    
Dx = eliminateBoundary(Dx, 1);
Dy = eliminateBoundary(Dy, 1);


function plotFlow(Dx, Dy, iter, iNum)
    cla;
    [yNum xNum] = size(Dx);
    [Y X]   = ndgrid(1:yNum, 1:xNum);
    sample  = ceil(yNum/45);
    IndexY  = 1:sample:yNum;
    IndexX  = 1:sample:xNum;
    scale   = sample*2;
    quiver(X(IndexY,IndexX),Y(IndexY,IndexX),...
           scale*Dx(IndexY,IndexX),scale*Dy(IndexY,IndexX),0,'-k');
    title(sprintf('Iteration %d of %d.',iter,iNum));
    axis ij equal; axis([-10 xNum+10 -10 yNum+10]);
    drawnow;
    