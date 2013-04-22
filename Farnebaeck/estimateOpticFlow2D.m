function [Dx Dy D] = estimateOpticFlow2D(ImgSeq, opt)
% estimateOpticFlow2D
%   ImgSeq  - Image sequence as a cube with dimensions: 
%             height x width x frames.
%   opt     - Struture with options:
%             * basisLen    - Length of Gaussian for basis vectors.
%             * basisSigma  - Standard deviation of Gaussian for basis.
%             * gamma       - Weight parameter between even and odd part of
%                             tensor.
%             * tensorLen   - Length of Gaussian for tensor smoothing.
%             * tensorSigma - Standard deviation of Gaussian for tensor.
%             * Confidence  - Matrix with confidence values which matches
%                             the dimensions of a frame of the image
%                             sequence.
%                             
%
% RETURN
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             This method computes the optic flow for the center frame of 
%             the image sequence, thus, Dx and Dy have the dimensions:
%             height x width.
%   D       - Distance between measurement and model estimate.
%
% DESCRIPTION
%   An implementation of 
%   Farnebaeck G. (2000). Fast and accurate motion estimation using 
%       orientation tensors and parametric motion models. In Proceedings 
%       of the 15th International Conference on Pattern Recognition 1, 
%       135-139.
%
%   Copyright (C) 2013  Florian Raudies, 01/04/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Get size parameters.
yNum = size(ImgSeq,1);
xNum = size(ImgSeq,2);
tNum = size(ImgSeq,3);
% Set default values for parameters of the method.
if nargin<2,                    opt             = struct(); end % UNITS
if ~isfield(opt,'basisLen'),    opt.basisLen    = 11;       end % pixels
if ~isfield(opt,'basisSigma'),  opt.basisSigma  = 1.6;      end % pixels
if ~isfield(opt,'gamma'),       opt.gamma       = 1/256;    end % -
if ~isfield(opt,'tensorLen'),   opt.tensorLen   = 41;       end % pixels
if ~isfield(opt,'tensorSigma'), opt.tensorSigma = 6.5;      end % pixels
if ~isfield(opt,'Confidence'),
    border = 5;
    opt.Confidence = ones(yNum,xNum);
    opt.Confidence(1:border,:) = 0;   opt.Confidence(end-border:end,:) = 0;
    opt.Confidence(:,1:border) = 0;   opt.Confidence(:,end-border:end) = 0;
end
% Retrieve parameter values.
basisLen    = opt.basisLen;
basisSigma  = opt.basisSigma;
gamma       = opt.gamma;
tensorLen   = opt.tensorLen;
tensorSigma = opt.tensorSigma;
Confidence  = opt.Confidence;
% Define weight matrix for tensor.
TensorW     = fspecial('gaussian', [tensorLen 1], tensorSigma);
basisLen2   = (basisLen-1)/2;
% Check if the provided sequence contains enough frames.
if tNum<basisLen, 
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], basisLen, tNum);
end
% Select only the center frame of the sequence within lsqConv3D.
opt.SelT = [ceil(tNum/2) ceil(tNum/2)] - basisLen2;
% Define weights.
A = exp(-(-basisLen2:+basisLen2).^2/(2*basisSigma^2));
% Compute minimum residual signal and motion tensor.
T = motionTensor(lsqConv3D(ImgSeq, A, opt), gamma);
% Remove the isotropic part from the motion tensor. Assumes tNum = 1.
T = removeIsotropic(squeeze(T));
% Define the constraint matrix for the affine motion model with six
% parameters. The constraint matrix is symmetric and we define each unique 
% entry once.
Q           = zeros(yNum, xNum, 25);
[Y X]       = ndgrid( 1:yNum, 1:xNum);
% Diagonal entries.
Q(:,:,1)    = X.*X.*T(:,:,1,1);
Q(:,:,2)    = Y.*Y.*T(:,:,1,1);
Q(:,:,3)    = T(:,:,1,1);
Q(:,:,4)    = X.*X.*T(:,:,2,2);
Q(:,:,5)    = Y.*Y.*T(:,:,2,2);
Q(:,:,6)    = T(:,:,2,2);
Q(:,:,7)    = T(:,:,3,3);
% Additional first row entries.
Q(:,:,8)    = X.*Y.*T(:,:,1,1);
Q(:,:,9)    = X.*T(:,:,1,1);
Q(:,:,10)   = X.*X.*T(:,:,1,2);
Q(:,:,11)   = X.*Y.*T(:,:,1,2);
Q(:,:,12)   = X.*T(:,:,1,2);
Q(:,:,13)   = X.*T(:,:,1,3);
% Additional second row entries.
Q(:,:,14)   = Y.*T(:,:,1,1);
Q(:,:,15)   = Y.*Y.*T(:,:,1,2);
Q(:,:,16)   = Y.*T(:,:,1,2);
Q(:,:,17)   = Y.*T(:,:,1,3);
% Additional third row entries (no symmetries and multiples).
Q(:,:,18)   = T(:,:,1,2);
Q(:,:,19)   = T(:,:,1,3);
% Additional fourth row entries.
Q(:,:,20)   = X.*Y.*T(:,:,2,2);
Q(:,:,21)   = X.*T(:,:,2,2);
Q(:,:,22)   = X.*T(:,:,2,3);
% Additional fifth row entries.
Q(:,:,23)   = Y.*T(:,:,2,2);
Q(:,:,24)   = Y.*T(:,:,2,3);
% Additional sixth row entries.
Q(:,:,25)   = T(:,:,2,3);
% Multpliciation with confidence values and normalization with confidence
% values in the neighborhood defined by TensorW.
Q = imfilter(imfilter(Q.*repmat(Confidence,[1 1 25]), TensorW), TensorW') ...
    ./( eps + repmat(imfilter(imfilter(Confidence, TensorW), TensorW'), [1 1 25]) );
% Define the full 7 x 7 matrix Q with all its entries for all yNum x xNum 
% positions.
Q = reshape(Q(:,:,[ 1  8  9 10 11 12 13;
                    8  2 14 11 15 16 17;
                    9 14  3 12 16 18 19;
                    10 11 12  4 20 21 22;
                    11 15 16 20  5 23 24;
                    12 16 18 21 23  6 25;
                    13 17 19 22 24 25  7]),[yNum xNum 7 7]);
% Solve the reduced system for each spatial location, see Eq. (14) in 
% Farnebaeck (2000).
J = Q(:,:,1:6,1:6);
B = Q(:,:,1:6,7);
P = zeros(yNum,xNum,6);
for ix=1:xNum,
    for iy=1:yNum,
        P(iy,ix,:) = - squeeze(J(iy,ix,:,:)) \ squeeze(B(iy,ix,:));
    end
end
% Compute the motion according to the estimated six parameters.
Dx = squeeze(sum(P(:,:,1:3).*cat(3, X,Y,ones(yNum,xNum)), 3));
Dy = squeeze(sum(P(:,:,4:6).*cat(3, X,Y,ones(yNum,xNum)), 3));
% Comput the minimum value, if requested.
if nargout > 2
    D = Q(:,:,7,7) + sum(P.*Q(:,:,1:6,7),3);
end
