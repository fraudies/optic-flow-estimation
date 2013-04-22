function [ImgSeq maxSpeed] = randomTextureExpansion(opt)
% randomTextureExpansion
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
%   Copyright (C) 2013  Florian Raudies, 01/02/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for parameters
if nargin<1,                opt         = struct(); end % UNITS
if ~isfield(opt,'h'),       opt.h       = 128;      end % pixels
if ~isfield(opt,'w'),       opt.w       = 128;      end % pixels
if ~isfield(opt,'fNum'),    opt.fNum    = 15;       end % frames
if ~isfield(opt,'expRate'), opt.expRate = 1.005;    end % units per frame
% Retrieve parameter valus.
h       = opt.h;
w       = opt.w;
fNum    = opt.fNum;
expRate = opt.expRate;
% Define random texture and image coordinates.
Texture = rand(h,w);
[Y X]   = ndgrid(linspace(-1,1,h), linspace(-1,1,w));
% Define image sequence and initalize.
ImgSeq  = zeros(h,w,fNum);
ImgSeq(:,:,1) = Texture;
% Add all other frames using an expansion of expRate and cubic
% interpolation.
for iFrame = 2:fNum,
    rate = expRate^(iFrame-1);
    ImgSeq(:,:,iFrame) = interp2(rate*X,rate*Y,Texture, X,Y, 'cubic');
end
% Calculate the maximum pixel speed as the maximum difference which is
% expRate-1 scaled to pixels by max(h,w)/2, since 1,...,w and 1,...,h map
% to the intervals -1...1 and -1...1, respectively. Thus, the maximum speed
% is 1 + max(h,w)/2*(expRate-1).
maxSpeed = 1 + max(h,w)/2*(expRate-1);
