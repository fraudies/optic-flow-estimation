function [ImgSeq maxSpeed] = randomTextures2linearMotions(opt)
% randomTextures2linearMotions
%   opt   - Structure with the following fields:
%           * h       - Height of video in pixels.
%           * w       - Width of video in pixels.
%           * fNum    - Number of frames of video.
%           * Vel1    - First velocity in pixels per frame.
%           * Vel2    - Second velocity in pixels per frame.
%
% RETURN
%   ImgSeq   - Image sequence as 3D matrix with dimension: h x w x fNum.
%   maxSpeed - Maximum speed in pixels per frame.
%
%   Copyright (C) 2013  Florian Raudies, 01/07/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Set default parameter values.
if nargin<1,                opt         = struct(); end % UNITS
if ~isfield(opt,'h'),       opt.h       = 128;      end % pixels
if ~isfield(opt,'w'),       opt.w       = 128;      end % pixels
if ~isfield(opt,'fNum'),    opt.fNum    = 15;       end % frames
if ~isfield(opt,'Vel1'),    opt.Vel1    = [0 -1];   end % pixel per frame
if ~isfield(opt,'Vel2'),    opt.Vel2t   = [0 +1];   end % pixel per frame
% Retrieve parameter values.
h       = opt.h;
w       = opt.w;
fNum    = opt.fNum;
Vel1    = opt.Vel1;
Vel2    = opt.Vel2;
% Define random texture and apply a Gaussian to reduce high frequencies.
Texture1 = imfilter(rand(h,w),...
            fspecial('gaussian',[5 5],1),'same','replicate');
Texture2 = imfilter(rand(h,w),...
            fspecial('gaussian',[5 5],1),'same','replicate');
% Define image sequence and initalize.
ImgSeq   = zeros(h,w,fNum);
% Add all other frames using cyclic boundary conditions.
for iFrame = 1:fNum,
    Delta1 = Vel1*(iFrame-1);
    Delta2 = Vel2*(iFrame-1);
    ImgSeq(:,:,iFrame) = 0.5 * (circshift(Texture1, Delta1) ...
                              + circshift(Texture2, Delta2));
end
maxSpeed = max(hypot(Vel1(1),Vel1(2)), hypot(Vel2(1),Vel2(2)));
