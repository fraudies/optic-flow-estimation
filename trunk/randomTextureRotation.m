function [ImgSeq maxSpeed] = randomTextureRotation(opt)
% randomTextureRotation
%   opt   - Structure with the following fields:
%           * h       - Height of video in pixels.
%           * w       - Width of video in pixels.
%           * fNum    - Number of frames of video.
%           * omega   - Rotation speed in radians per frame.
%
% RETURN
%   ImgSeq   - Image sequence as 3D matrix with dimension: h x w x fNum.
%   maxSpeed - Maximum speed in pixels per frame.
%
%   Copyright (C) 2013  Florian Raudies, 01/07/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default values for parameters
if nargin<1,                opt         = struct(); end % UNITS
if ~isfield(opt,'h'),       opt.h       = 128;      end % pixels
if ~isfield(opt,'w'),       opt.w       = 128;      end % pixels
if ~isfield(opt,'fNum'),    opt.fNum    = 15;       end % frames
if ~isfield(opt,'omega'),   opt.omega   = pi/15;    end % radians per frame
% Retrieve prameters.
h       = opt.h;
w       = opt.w;
fNum    = opt.fNum;
omega   = opt.omega;
rate    = sqrt(2);
oH      = ceil(h*rate);
oH      = oH + mod(oH,2);
oW      = ceil(w*rate);
oW      = oW + mod(oW,2);
SelY    = ceil(oH/2) + (1:h) -h/2;
SelX    = ceil(oW/2) + (1:w) -w/2;
% Define random texture and coordinates.
Texture = rand(oH,oW);
[Y X]   = ndgrid(linspace(-1,1,oH), linspace(-1,1,oW));
% Define image sequence and initalize.
ImgSeq  = zeros(h,w,fNum);
% Add all other frames using an expansion of expRate and cubic
% interpolation.
for iFrame = 1:fNum,
    theta = omega * (iFrame - 1);
    Xt = X*cos(theta)-Y*sin(theta);
    Yt = X*sin(theta)+Y*cos(theta);
    Img = interp2(X,Y, Texture, Xt,Yt, 'cubic');
    ImgSeq(:,:,iFrame) = Img(SelY,SelX);
end
% Calculate the maximum pixel speed as the maximum difference which is
% 1*omega which is scaled to pixels by max(h,w)/2, since 1,...,w and 
% 1,...,h map to the intervals -1...1 and -1...1, respectively. Thus, the 
% maximum speed is 1 + max(h,w)/2*omega.
maxSpeed = 1 + max(h,w)/2*omega;
ImgSeq(isnan(ImgSeq)) = 0;
