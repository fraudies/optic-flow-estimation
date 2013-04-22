function [Dx Dy] = estimateOpticFlow2D(ImgSeq, opt)
% estimateOpticFlow2D
%   ImgSeq  - Image sequence as a cube with dimensions: 
%             height x width x frames.
%   opt     - Struture with options:
%             * k       - Parameter of temporal filter kernel.
%             * htNum   - Number of samples for temporal filter kernel.
%             * n1,n2   - Parameter(s) of temporal filter kernel for early
%                         and futher delayed response.
%             * oNum    - Number of orientations for spatial Gabors.
%             * f       - Spatial frequency of Gabor filter in cycles per 
%                         pixel.
%             * sigma   - Standard deviation of Gabor filter.
%
% RETURN
%   Dx      - X-component of optic flow.
%   Dy      - Y-component of optic flow.
%             Computes flow for an interpolated frame pair that is
%             determined by the parameters k, n1, and n2 for the temporal
%             filter kernel. The computed flow Dx and Dy has the
%             dimensions: height x width.
%
% DESCRIPTION
%   An implementation of the idea of 
%   Adelson, E.H. and Bergen, J.R. (1985). Spatiotemporal energy models for 
%       the perception of motion. Journal of the Optical Society of 
%       America 2(2), 284-299.
%   I did not add a sampling of different speeds, which could be realized
%   e.g. by changing the parameter k. I use the orientation of the Gabor
%   filter as indication of the motion direction and, thus, detect only
%   motion normal to the Gabor filter's orientation.
%
%   Copyright (C) 2013  Florian Raudies, 01/08/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for paraemters of the method.
if nargin<2,                opt         = struct(); end
if ~isfield(opt,'k'),       opt.k       = 0.1;      end
if ~isfield(opt,'htNum'),   opt.htNum   = 15;       end
if ~isfield(opt,'n1'),      opt.n1      = 3;        end
if ~isfield(opt,'n2'),      opt.n2      = 5;        end
if ~isfield(opt,'oNum'),    opt.oNum    = 8;        end
if ~isfield(opt,'f'),       opt.f       = 1/4;      end % cycles per pixel
if ~isfield(opt,'sigma'),   opt.sigma   = 4;        end % pixels
% Retrieve parameter values.
k       = opt.k;
htNum   = opt.htNum;
oNum    = opt.oNum;
n1      = opt.n1;
n2      = opt.n2;
sigma   = opt.sigma;
f       = opt.f;
% Check if the provided sequence contains enough frames.
tNum    = size(ImgSeq,3);
if tNum<htNum,
    error('MATLAB:frameErr', ['This method requires at least %d frames ',...
        'but only %d frames were provided!'], htNum, tNum);
end
% Define causal temporal filter kernel.
hi      = @(T,k,ni) (k*T).^ni.*exp(-k*T.^2) ...
                  .*(1/factorial(ni)-(k*T).^2/factorial(ni+2));
T       = 1:htNum;
H1      = hi(T,k,n1);
H2      = hi(T,k,n2);
% Normalize filter kernels.
H1      = H1/sum(H1);
H2      = H2/sum(H2);
Ori     = (1:oNum)/oNum*pi;
yNum    = size(ImgSeq,1);
xNum    = size(ImgSeq,2);
fKernel     = @(sigma)      (-floor(3*sigma):+1:+floor(3*sigma));
fGaussian   = @(Z,sigma)    1/(sqrt(2*pi)*sigma) * exp(-Z.^2/(2*sigma^2));
fGaborSin   = @(Z,f,sigma)  fGaussian(Z,sigma) .* sin(2*pi*f*Z);
fGaborCos   = @(Z,f,sigma)  fGaussian(Z,sigma) .* cos(2*pi*f*Z);
Kernel      = fKernel(sigma);
kNum        = length(Kernel);
SelT     = 1:htNum;
ImgSeq   = ImgSeq(:,:,SelT);
% Causal filtering in the temporal domain.
ImgSeqH1 = squeeze(sum(ImgSeq.*shiftdim(repmat(H1(:),[1 yNum xNum]),1),3));
ImgSeqH2 = squeeze(sum(ImgSeq.*shiftdim(repmat(H2(:),[1 yNum xNum]),1),3));
R1 = zeros(yNum,xNum,oNum);
R2 = zeros(yNum,xNum,oNum);
for iOri = 1:oNum,
    fx = f*cos(Ori(iOri));
    fy = f*sin(Ori(iOri));
    Gabor.SinX = reshape(fGaborSin(Kernel,fx,sigma), [1 kNum]);
    Gabor.CosX = reshape(fGaborCos(Kernel,fx,sigma), [1 kNum]);
    Gabor.SinY = reshape(fGaborSin(Kernel,fy,sigma), [kNum 1]);
    Gabor.CosY = reshape(fGaborCos(Kernel,fy,sigma), [kNum 1]);
    R1(:,:,iOri) = filterSepGabor(ImgSeqH1, Gabor);
    R2(:,:,iOri) = filterSepGabor(ImgSeqH2, Gabor);
end
% Normalization over orientations.
R1 = R1./(eps + repmat(mean(abs(R1),3),[1 1 oNum]));
R2 = R2./(eps + repmat(mean(abs(R2),3),[1 1 oNum]));
% Calculate opponent motion energy.
E = 4*(real(R1).*imag(R2) - real(R2).*imag(R1));
% Negative sign encodes the opposite direction.
E = max(cat(3, -E, E),0);
% Define prototopye motion vectors Vx and Vy.
Vx = cos([Ori pi+Ori]);
Vy = sin([Ori pi+Ori]);
% Return vector corresponding to maximum motion energy.
[~, Index] = max(E,[],3);
Dx = Vx(Index);
Dy = Vy(Index);

function ImgSeq = filterSepGabor(ImgSeq, Gabor)
% Calculate filter response for 2D spatial Gabors.
GaborCosX = imfilter(ImgSeq, Gabor.CosX, 'same', 'replicate');
GaborSinX = imfilter(ImgSeq, Gabor.SinX, 'same', 'replicate');
ImgSeq = imfilter(GaborSinX, Gabor.CosY, 'same', 'replicate') ...
       + imfilter(GaborCosX, Gabor.SinY, 'same', 'replicate') ...
       + 1i ...
       *(imfilter(GaborCosX, Gabor.CosY, 'same', 'replicate') ...
       - imfilter(GaborSinX, Gabor.SinY, 'same', 'replicate'));    