function [ImgSeq maxSpeed] = randomTexturesExpCon(opt)
% randomTexturesExpCon
%   opt   - Structure with the following fields:
%           * h       - Height of video in pixels.
%           * w       - Width of video in pixels.
%           * fNum    - Number of frames of video.
%           * expRate - Expansion rate. Choose a small rate, e.g. 1.005.
%
% RETURN
%   ImgSeq   - Image sequence as 3D matrix with dimension: h x w x fNum.
%   maxSpeed - Maximum speed in pixels per frame.
%
%   Copyright (C) 2013  Florian Raudies, 01/07/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default parameters values.
if nargin<1,                opt         = struct(); end % UNITS
if ~isfield(opt,'h'),       opt.h       = 128;      end % pixels
if ~isfield(opt,'w'),       opt.w       = 128;      end % pixels
if ~isfield(opt,'fNum'),    opt.fNum    = 15;       end % frames
if ~isfield(opt,'expRate'), opt.expRate = 1.005;    end % units per frame
% Retrieve parameter values.
h       = opt.h;
w       = opt.w;
fNum    = opt.fNum;
expRate = opt.expRate;
rate    = expRate^(fNum-1);
oH      = ceil(h*rate);
oH      = oH + mod(oH,2);
oW      = ceil(w*rate);
oW      = oW + mod(oW,2);
SelY    = ceil(oH/2) + (1:h) -h/2;
SelX    = ceil(oW/2) + (1:w) -w/2;
% Define random texture and apply a Gaussian to reduce high frequencies.
Texture1 = imfilter(rand(oH,oW),...
            fspecial('gaussian',[5 5],1),'same','replicate');
Texture2 = imfilter(rand(oH,oW),...
            fspecial('gaussian',[5 5],1),'same','replicate');
[Y X]   = ndgrid(linspace(-1,1,oH), linspace(-1,1,oW));
% Define image sequence and initalize.
ImgSeq  = zeros(h,w,fNum);
% Add all other frames using an expansion of expRate and cubic
% interpolation.
for iFrame = 1:fNum,
    rate = expRate^(iFrame-1);
    Img1 = interp2(rate*X,rate*Y,Texture1, X,Y, 'cubic');
    Img2 = interp2(X,Y,Texture2, rate*X,rate*Y, 'cubic');
    ImgSeq(:,:,iFrame) = 0.5 * (Img1(SelY,SelX) + Img2(SelY,SelX));
end
% Calculate the maximum pixel speed as the maximum difference which is
% expRate-1 scaled to pixels by max(h,w)/2, since 1,...,w and 1,...,h map
% to the intervals -1...1 and -1...1, respectively. Thus, the maximum speed
% is 1 + max(h,w)/2*(expRate-1).
maxSpeed = 1 + max(h,w)/2*(expRate-1);
