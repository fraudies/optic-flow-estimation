function ImgSeq = readImgSeq(filePattern, firstFrame, lastFrame)
% readImgSeq
%   filePattern - This describes file path and name pattern to be used in 
%                 the sprintf command together with the current frame 
%                 number to read each frame of the image sequence.
%   firstFrame  - Index for the first frame.
%   lastFrame   - Index for the last frame.
%
% RETURN
%   ImgSeq      - 3D image cube with dimensions: height x width x frames.
%
% DESCRITPION
%   The method reads (lastFrame-firstFrame+1) image frames and stores them 
%   in a 3D image cube.
% 
%   Copyright (C) 2013  Florian Raudies, 01/14/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.

fNum        = lastFrame - firstFrame + 1; % number of frames
[yNum xNum] = size(imread(sprintf(filePattern,0))); % get dimensions
ImgSeq      = zeros(yNum,xNum,fNum);
for iFrame = 1:fNum,
    ImgSeq(:,:,iFrame) = double(imread(...
        sprintf(filePattern, firstFrame + iFrame - 1)))/255;
end
