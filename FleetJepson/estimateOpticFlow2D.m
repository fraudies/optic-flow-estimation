function [Dx Dy] = estimateOpticFlow2D(ImgSeq, opt)
% estimateOpticFlow2D
%   ImgSeq  - Image sequence as a cube with dimensions: 
%             height x width x frames.
%   opt     - Struture with options:
%             * f       - Spatio-temporal frequency in cycles per pixel or
%                         cycles per frame. The spatial frequency f is the 
%                         radius of ceil(pi/2*(2^beta+1)/(2^beta-1)) tested
%                         orientations.
%             * beta    - Frequency bandwidth in octaves.
%             * FreqTmp - Temporal frequencies.
%             In total a combination of ceil(pi/2*(2^beta+1)/(2^beta-1))
%             times length(FreqTmp) filter responses are calculated and
%             evaluated.
%
% RETURN
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             Computes flow for the center frame of the sequence, thus, Dx
%             and Dy are of dimensions: height x width.
%
% DESCRIPTION
%   An implementation of the idea of 
%   Fleet, D. and Jepson, A.D. (1990). Computation of component image
%       velocity from local phase information. International Journal of 
%       Computer Vision 5(1), 77-104.
%   Note that the paper describes most parameters of the algorithm on p. 86
%   as speed values, number of directions, and directional increment.
%   However, from the paper it is unlcear how these speed values link to
%   Gabor filters. For instance, using ds = pi/sigma is not defined for
%   zero speeds.
%   A shorter and sometimes clearer write-up of the method is given in
%   Fleet, D.J. and Jepson, A.D. (1989). Computation of normal velocity 
%       from local phase information. In Proceedings of Computer Vision and 
%       Pattern Recognition, CVPR'89. 379-386.
%   These papers do not describe how to combine responses from different 
%   frequencies and in particular the normalization accross them. Here, I 
%   weigh estimates by their amplitude, which is later re-normlized to get 
%   the correct speed (see lines 108 and 131 & 132).
%
%   Copyright (C) 2013  Florian Raudies, 05/16/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for paraemters of the method.
if nargin<2,                opt         = struct();         end % UNITS
if ~isfield(opt,'f'),       opt.f       = 1/4;              end % cycles per pixel or frame
if ~isfield(opt,'beta'),    opt.beta    = 0.8;              end % octaves
if ~isfield(opt,'FreqTmp'), opt.FreqTmp = [-opt.f 0 opt.f]; end % temporal frequencies
% Retrieve parameters.
f           = opt.f;
beta        = opt.beta;
FreqTmp     = opt.FreqTmp;
% Define additional prameters and Gabor functions.
sigmaFreq   = f*(2^beta-1)/(2^beta+1);  % sigma in Fourier domain
sigma       = 1/(2*pi*sigmaFreq);       % sigma in non-transformed domain
nOri        = ceil(pi/2*(2^beta+1)/(2^beta-1)); % number of orientations
Ori         = (1:nOri)/nOri*pi;
b           = (2^beta+1)/(2^beta-1);    % parameter to remove DC part of even Gabor
fKernel     = @(sigma)      (-floor(3*sigma):+1:+floor(3*sigma));
fGaussian   = @(Z,sigma)    1/(sqrt(2*pi)*sigma) * exp(-Z.^2/(2*sigma^2));
fGaborSin   = @(Z,f,sigma)  fGaussian(Z,sigma) .* sin(2*pi*f*Z);
% DC part removed. Note, the cosine/sine relationship is not effected much.
fGaborCos   = @(Z,f,sigma,b)fGaussian(Z,sigma) .* (cos(2*pi*f*Z)-exp(-b^2/2));
Kernel      = fKernel(sigma);
kNum        = length(Kernel);
yNum        = size(ImgSeq,1);
xNum        = size(ImgSeq,2);
nTmp        = length(FreqTmp);  % number of temporal frequencies
D           = [-1 0 1]/2;       % differential operator
dNum        = length(D);
% Define the space for angular frequencies in x, y, and t. These are named 
% kx, ky, and omega, respectively.
Kx      = zeros(yNum,xNum,nOri,nTmp); % spatial frequency
Ky      = zeros(yNum,xNum,nOri,nTmp); % spatial frequency
Omega   = zeros(yNum,xNum,nOri,nTmp); % temporal frequency
A       = zeros(yNum,xNum,nOri,nTmp); % amplitude signal
% Check if the provided sequence contains enough frames.
tNum    = size(ImgSeq,3);
if tNum<(kNum+2), 
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], kNum+2, tNum);
end
% Select the center number of kNum+2 frames.
center  = ceil(tNum/2);
SelT    = center + ( -(kNum-1)/2-1 : +(kNum-1)/2+1 ); % +-1 for D.
ImgSeq  = ImgSeq(:,:,SelT);
% Compute the angular frequencies in spatial and temporal domain.
for iTmp = 1:nTmp,
    ft = FreqTmp(iTmp);
    Gabor.SinT = reshape(fGaborSin(Kernel,ft,sigma),   [1 1 kNum]);
    Gabor.CosT = reshape(fGaborCos(Kernel,ft,sigma,b), [1 1 kNum]);
    for iOri = 1:nOri,
        fx = f*cos(Ori(iOri));
        fy = f*sin(Ori(iOri));
        Gabor.SinX = reshape(fGaborSin(Kernel,fx,sigma),   [1 kNum 1]);
        Gabor.CosX = reshape(fGaborCos(Kernel,fx,sigma,b), [1 kNum 1]);
        Gabor.SinY = reshape(fGaborSin(Kernel,fy,sigma),   [kNum 1 1]);
        Gabor.CosY = reshape(fGaborCos(Kernel,fy,sigma,b), [kNum 1 1]);
        % Compute the Gabor filter response.
        R   = filterSepGabor(ImgSeq, Gabor);
        % Compute the temoral derivative.
        Rt  = convn(R, reshape(D, [1 1 dNum]), 'valid');
        % Select the center frame.
        R   = squeeze(R(:,:, 2));
        % Compute the spatial derivatives.
        Rx  = imfilter(R, reshape(D, [1 dNum 1]), 'same', 'replicate');
        Ry  = imfilter(R, reshape(D, [dNum 1 1]), 'same', 'replicate');
        RInv = 1./(abs(R) + eps);
        % Compute the spatial and temporal phases and add them to the 
        % tested ones. See Fleet & Jepson (1989), Eq. (6) on page 380.
        Kx(:,:,iOri,iTmp)    = (imag(Rx).*real(R) ...
                               -real(Rx).*imag(R)).*RInv + 2*pi*fx;
        Ky(:,:,iOri,iTmp)    = (imag(Ry).*real(R) ...
                               -real(Ry).*imag(R)).*RInv + 2*pi*fy;
        Omega(:,:,iOri,iTmp) = (imag(Rt).*real(R) ...
                               -real(Rt).*imag(R)).*RInv + 2*pi*ft;
        A(:,:,iOri,iTmp)     = abs(R);
    end
end
% Put all constraints into a single dimension.
Kx      = reshape(Kx,    [yNum xNum nOri*nTmp]);
Ky      = reshape(Ky,    [yNum xNum nOri*nTmp]);
Omega   = reshape(Omega, [yNum xNum nOri*nTmp]);
A       = reshape(A,     [yNum xNum nOri*nTmp]);
% Compute the normal vector.
LenInv  = 1./(hypot(Kx, Ky) + eps);
Nx      = +Kx.*LenInv;      % x-component of normal
Ny      = +Ky.*LenInv;      % y-component of normal
S       = -Omega.*LenInv;   % speed along normal
% Compute normal flow. Noramlization using the amplitude.
Dx      = squeeze(mean(S.*Nx,3)./mean(A,3));
Dy      = squeeze(mean(S.*Ny,3)./mean(A,3));


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
