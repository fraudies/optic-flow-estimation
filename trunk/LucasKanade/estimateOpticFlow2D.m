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
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             Computes the optic flow for all frames of the sequence,
%             thus, Dx and Dy have dimensions: height x width x frames-1.
%
% DESCRIPTION
%   A modern implementation of the idea of 
%   Lucas, B.D. and Kanade, T. (1981). An iterative image registration 
%       technique with and application to stereo vision. In Proceedings 
%       of Imaging Understanding Workshop, 121-130.
%
%   Copyright (C) 2013  Florian Raudies, 12/31/2012, Boston University.
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
% Compute the partial derivatives in x, y, and t.
ktNum    = size(DiffT,3);
ImgSeqDx = imfilter(ImgSeq, DiffX, 'same', 'replicate');
ImgSeqDy = imfilter(ImgSeq, DiffY, 'same', 'replicate');
ImgSeqDt = convn(ImgSeq, DiffT, 'valid');
% Select the spatial and temporal valid part of the partial derivatives.
ValidT   = (1+floor((ktNum-1)/2)) : (tNum-floor(ktNum/2));
ImgSeqDx = squeeze(ImgSeqDx(:,:,ValidT));
ImgSeqDy = squeeze(ImgSeqDy(:,:,ValidT));
% Compute the coefficients A, B, and C.
A = imfilter(ImgSeqDx.^2,       W, 'same', 'replicate');
B = imfilter(ImgSeqDx.*ImgSeqDy,W, 'same', 'replicate');
C = imfilter(ImgSeqDy.^2,       W, 'same', 'replicate');
D = -1./(eps+A.*C-B.^2);
E = imfilter(ImgSeqDx.*ImgSeqDt,W, 'same', 'replicate');
F = imfilter(ImgSeqDy.*ImgSeqDt,W, 'same', 'replicate');
% Compute the flow variables.
Dx = D.*(+A.*E -B.*F);
Dy = D.*(-B.*E +C.*F);
