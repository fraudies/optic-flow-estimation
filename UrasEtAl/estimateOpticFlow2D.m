function [Dx Dy Valid] = estimateOpticFlow2D(ImgSeq, opt)
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
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             Computes the optic flow for all frames of the sequence,
%             thus, Dx and Dy have dimensions: height x width x frames-1.
%
% DESCRIPTION
%   An implementation of the idea of 
%   Uras, S., Girosi, F., Verri, A. and Torre, V. (1988). A computational 
%       approach to motion perception. Biological Cybernetics 60, 79-87.
%   I did not apply any post-processing to the optic flow field as e.g.
%   described on page 81, right column. I also use a different method for 
%   computing the validity of the estimated flow.
%
%   Copyright (C) 2013  Florian Raudies, 01/04/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for paraemters of the method.
if nargin<2,                opt         = struct(); end
if ~isfield(opt,'W'),       opt.W       = fspecial('gaussian',[9 9],2); end
if ~isfield(opt,'DiffX'),   opt.DiffX   = [-1 8 0 -8 1]/12; end
if ~isfield(opt,'DiffY'),   opt.DiffY   = opt.DiffX'; end
if ~isfield(opt,'DiffT'),   opt.DiffT   = reshape([-1 1],[1 1 2]); end
% Retreive default parmeters to be saved in their own variables.
W       = opt.W;
DiffX   = opt.DiffX;
DiffY   = opt.DiffY;
DiffT   = opt.DiffT;
% Check if the provided sequence contains at least two frames.
tNum     = size(ImgSeq,3);
if tNum<2, 
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], 2, tNum);
end
% Compute the Hessian matrix H = [d^2f/dx^2 d^2f/dxdy]  = [H11 H12]
%                                [d^2f/dxdy d^2f/dy^2].   [H12 H22]
tNum     = size(ImgSeq,3);
ktNum    = size(DiffT,3);
ImgSeqDx = imfilter(ImgSeq, DiffX, 'same', 'replicate');
ImgSeqDy = imfilter(ImgSeq, DiffY, 'same', 'replicate');
ImgSeqDt = convn(ImgSeq, DiffT, 'valid');
% Select the spatial and temporal valid part of the partial derivatives.
ValidT   = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
ImgSeqDx = squeeze(ImgSeqDx(:,:,ValidT));
ImgSeqDy = squeeze(ImgSeqDy(:,:,ValidT));
% Hessian matrix.
H11 = imfilter(ImgSeqDx, DiffX, 'same', 'replicate');
H12 = imfilter(ImgSeqDx, DiffY, 'same', 'replicate'); % = H21
H22 = imfilter(ImgSeqDy, DiffY, 'same', 'replicate');
% Right-hand side.
B1  = imfilter(ImgSeqDt, DiffX, 'same', 'replicate');
B2  = imfilter(ImgSeqDt, DiffY, 'same', 'replicate');
% Compute local neighborhood.
H11 = imfilter(H11, W, 'same', 'replicate');
H12 = imfilter(H12, W, 'same', 'replicate');
H22 = imfilter(H22, W, 'same', 'replicate');
B1  = imfilter(B1,  W, 'same', 'replicate');
B2  = imfilter(B2,  W, 'same', 'replicate');
% Compute the flow variables.
D  = -1./(H11.*H22 - H12.^2 + eps);
Dx = D.*(+H22.*B1 -H12.*B2);
Dy = D.*(-H12.*B1 +H11.*B2);
% Validation of estimates and assumption about of flow from a plane.
% Compute the matrix M = [du/dx du/dy] which is the derivative of flow.
%                        [dv/dx dv/dy]
M11 = imfilter(Dx, DiffX, 'same', 'replicate');
M12 = imfilter(Dx, DiffY, 'same', 'replicate');
M21 = imfilter(Dy, DiffX, 'same', 'replicate');
M22 = imfilter(Dy, DiffY, 'same', 'replicate');
% % It seems that the normalization with hypot(B1,B2) does not lead to the
% % desired criteria. Note that M*dE =! 0.
% Delta = hypot(M11.*ImgSeqDx + M21.*ImgSeqDy, ...
%               M12.*ImgSeqDx + M22.*ImgSeqDy); % ./ (hypot(B1, B2) + eps);
% Valid = Delta < 0.05; % Valid with respect to assumption of const. motion.
% An alternative criteria of numerical stability can be formulated using 
% the condition number of the Hessian matrix.
% % Define validity of solution through the condition number of the Hessian.
% DiscriminantSqrt = sqrt(4*H12.^2 + (H11-H22).^2);
% LambdaPlus  = 0.5 * (H11+H22) + DiscriminantSqrt;
% LambdaMinus = 0.5 * (H11+H22) - DiscriminantSqrt;
% CondOfH = LambdaPlus./(LambdaMinus + eps);
% Valid = CondOfH < 0.1;
% Yet another criteria, and the one that is most intuitive, is to compute
% the hypothetical flow, which was left out using the current solution and 
% applying a threshold to this hypthetical flow part.
D1 = M11.*ImgSeqDx + M21.*ImgSeqDy;
D2 = M12.*ImgSeqDx + M22.*ImgSeqDy;
Delta = hypot(D.*(+H22.*D1 -H12.*D2), D.*(-H12.*D1 +H11.*D2));
Valid = Delta < 10;

