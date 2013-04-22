function [Dx Dy C] = estimateOpticFlow2D(ImgSeq, opt)
% estimateOpticFlow2D
%   ImgSeq  - Image sequence as a cube with dimensions:
%             height x width x frames.
%   opt     - Struture with options:
%             * sigma   - Standard deviation of Gaussian for spatial and
%                         temporal envelope.
%             * fMax    - Maximum frequency of sampled interval and -fMax
%                         is the minimum frequency.
%             * fSample - Sampling period in the frequency domain.
%
% RETURN
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             Both flow components of the dimensions: height x width x 2.
%   C       - Confidence values.
%
% DESCRIPTION
%   An implementation of
%   Shizawa, M. and Mase, K. (1991). A unified computational theory for
%       motion transparency and motion boundaries based on eigenergy
%       analysis. In Proceedings of IEEE Computer Vision and Pattern
%       Recognition Conference CVPR'91, 289-295.
%   I assume that two motion component are present in the input. Unlike the
%   original implementation, I use a different range of frequencies as
%   default. In addition I separately perform Gabor filtering in each 
%   domain and compute motions for each pixel.
%
%   Copyright (C) 2013  Florian Raudies, 01/07/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for parameters of the method.
if nargin<2,                opt         = struct(); end % UNITS
if ~isfield(opt,'sigma'),   opt.sigma   = 4;        end % pixels or frames
if ~isfield(opt,'fMax'),    opt.fMax    = 0.4;      end % cycles per pixel or frame
if ~isfield(opt,'fSample'), opt.fSample = 0.2;      end % cycles per pixel or frame
% Retrieve parameter values.
sigma   = opt.sigma;
fMax    = opt.fMax;
fSample = opt.fSample;
% Define spatio-temporal frequency space.
F           = -fMax : fSample : +fMax;
[Fy Fx Ft]  = ndgrid(F, F, F); % Define frequencies on a cubic lattice.
Fx          = Fx(:);
Fy          = Fy(:);
Ft          = Ft(:);
fNum        = numel(Fx);
% Define Gabor filters.
fKernel     = @(sigma)      (-floor(3*sigma):+1:+floor(3*sigma));
fGaussian   = @(Z,sigma)    1/(sqrt(2*pi)*sigma) * exp(-Z.^2/(2*sigma^2));
fGaborSin   = @(Z,f,sigma)  fGaussian(Z,sigma) .* sin(2*pi*f*Z);
fGaborCos   = @(Z,f,sigma)  fGaussian(Z,sigma) .* cos(2*pi*f*Z);
Kernel      = fKernel(sigma);
kNum        = length(Kernel);
yNum        = size(ImgSeq,1);
xNum        = size(ImgSeq,2);
% Mapping of indices: xx = 11, yy = 22, tt = 33, xy = 12, yt = 23, xt = 13.
SumASqFreq  = zeros(yNum,xNum,6,6);
SumASq      = zeros(yNum,xNum);
tNum        = size(ImgSeq,3);
% Check if the provided sequence contains enough frames.
if tNum<kNum, 
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], kNum, tNum);
end
% Trim the image sequence to a number of frames that leaves a single frame
% after all temporal filtering. This number is kNum.
center  = ceil(tNum/2);
SelT    = center + (-floor(kNum/2) : +1 : +floor(kNum/2));
ImgSeq  = ImgSeq(:,:,SelT);
for iFreq = 1:fNum,
    fx = Fx(iFreq);
    fy = Fy(iFreq);
    ft = Ft(iFreq);
    Gabor.SinT = reshape(fGaborSin(Kernel,ft,sigma), [1 1 kNum]);
    Gabor.CosT = reshape(fGaborCos(Kernel,ft,sigma), [1 1 kNum]);
    Gabor.SinX = reshape(fGaborSin(Kernel,fx,sigma), [1 kNum 1]);
    Gabor.CosX = reshape(fGaborCos(Kernel,fx,sigma), [1 kNum 1]);
    Gabor.SinY = reshape(fGaborSin(Kernel,fy,sigma), [kNum 1 1]);
    Gabor.CosY = reshape(fGaborCos(Kernel,fy,sigma), [kNum 1 1]);
    ASq     = abs(filterSepGabor(ImgSeq, Gabor)).^2;
    SumASq  = SumASq + ASq;
    FreqVec = [fx*fx fy*fy ft*ft fx*fy fy*ft fx*ft];
    for ii = 1:6,
        for jj = 1:6,
            SumASqFreq(:,:,ii,jj) = SumASqFreq(:,:,ii,jj) ...
                                  + ASq*FreqVec(ii)*FreqVec(jj);
        end
    end
end
% Define noramlized moments.
M = SumASqFreq./(eps + repmat(SumASq,[1 1 6 6]));
M(:,:,[4 5 6],:) = 2*M(:,:,[4 5 6],:);
M(:,:,:,[4 5 6]) = 2*M(:,:,:,[4 5 6]);
% Allocate space for velocities.
Dx  = zeros(yNum,xNum,2);
Dy  = zeros(yNum,xNum,2);
C   = zeros(yNum,xNum); % Confidence
% Compute the velocities.
for iy = 1:yNum,
    for ix = 1:xNum,
        Mxy = squeeze(M(iy,ix,:,:));
        [Vec Val] = eig(Mxy);
        [~, index] = min(Val(eye));
        Vs = Vec(:,index);
        V1V2 = [Vs(1) Vs(4) Vs(6); ...
                Vs(4) Vs(2) Vs(5); ...
                Vs(6) Vs(5) Vs(3)];
        [Vec Val] = eig(V1V2);
        % Val contains the sorted eigenvalues from smallest to largest.
        l1 = Val(1,1);  l2 = Val(2,2);  l3 = Val(3,3);
        Em = Vec(:,1);  Ep = Vec(:,3);
        lp = l1 - l2;
        lm = l3 - l2;
        V1 = +sqrt( -lm/(lp-lm) )*Em + sqrt( +lp/(lp-lm) )*Ep;
        V2 = -sqrt( -lm/(lp-lm) )*Em + sqrt( +lp/(lp-lm) )*Ep;
        % Multiply with sign of third velocity component to resolve the
        % ambiguity of the direction of the eigenvectors.
        Dx(iy,ix,[1 2]) = [V1(1)*sign(V1(3)) V2(1)*sign(V2(3))];
        Dy(iy,ix,[1 2]) = [V1(2)*sign(V1(3)) V2(2)*sign(V2(3))];
        C(iy,ix) = l2;
    end
end
C = abs(C);

function ImgSeq = filterSepGabor(ImgSeq, Gabor)
% Calculate 3D separable Gabor filter results.
GaborSinT = convn(ImgSeq, -Gabor.SinT, 'valid'); % '-' for convn
GaborCosT = convn(ImgSeq, +Gabor.CosT, 'valid');
GaborSinTCosX = imfilter(GaborSinT, Gabor.CosX, 'same', 'replicate');
GaborSinTSinX = imfilter(GaborSinT, Gabor.SinX, 'same', 'replicate');
GaborCosTSinX = imfilter(GaborCosT, Gabor.SinX, 'same', 'replicate');
GaborCosTCosX = imfilter(GaborCosT, Gabor.CosX, 'same', 'replicate');
ImgSeq  = imfilter(GaborCosTCosX, Gabor.CosY, 'same', 'replicate') ...
        - imfilter(GaborCosTSinX, Gabor.SinY, 'same', 'replicate') ...
        - imfilter(GaborSinTSinX, Gabor.CosY, 'same', 'replicate') ...
        - imfilter(GaborSinTCosX, Gabor.SinY, 'same', 'replicate') ...
        + 1i...
        *(imfilter(GaborSinTCosX, Gabor.CosY, 'same', 'replicate') ...
        - imfilter(GaborSinTSinX, Gabor.SinY, 'same', 'replicate') ...
        + imfilter(GaborCosTSinX, Gabor.CosY, 'same', 'replicate') ...
        + imfilter(GaborCosTCosX, Gabor.SinY, 'same', 'replicate'));