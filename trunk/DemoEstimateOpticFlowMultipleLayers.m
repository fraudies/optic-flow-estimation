clc
clear all
close all

% *************************************************************************
% Estimate multiple layers of motion (here two) for three image sequences.
% One is composed of two linera motions, another of expansion/contraction,
% and a third of clockwise and counterclockwise rotation.
%
% Attention: This script may take tens of seconds to finish!
%
%   Copyright (C) 2013  Florian Raudies, 05/16/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Select one of the two algorithms.
% algoPath = './MotaEtAl/';
algoPath = './ShizawaMase/';


% Generate the image sequences.
opt.h       = 128; % Set parameters for stimulus.
opt.w       = 128;
opt.fNum    = 27;
opt.Vel1    = [1 1];
opt.Vel2    = [0 -1];
[ImgSeq1 maxSpeed1] = randomTextures2linearMotions(opt);
opt.expRate = 1.02;
[ImgSeq2 maxSpeed2] = randomTexturesExpCon(opt);
opt.omega   = pi/100;
[ImgSeq3 maxSpeed3] = randomTextures2Rotations(opt);


% Display the estimated optic flow.
h       = size(ImgSeq1,1);
w       = size(ImgSeq1,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h; 
% For the display the flow is scaled by division with its maximum speed and
% multiplication with the sampling factor.
len1    = sample/maxSpeed1;
len2    = sample/maxSpeed2;
len3    = sample/maxSpeed3;


% *************************************************************************
% Estimate and display flow for the 1st sequence.
% *************************************************************************
addpath(algoPath); 
[Dx1 Dy1] = estimateOpticFlow2D(ImgSeq1);
rmpath(algoPath);
figure('Position',[50 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx1(IndexY,IndexX,1)*len1, Dy1(IndexY,IndexX,1)*len1,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx1(IndexY,IndexX,2)*len1, Dy1(IndexY,IndexX,2)*len1,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Linear motions.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;


% *************************************************************************
% Estimate and display flow for the 2nd sequence.
% *************************************************************************
addpath(algoPath); 
[Dx2 Dy2] = estimateOpticFlow2D(ImgSeq2);
rmpath(algoPath);
figure('Position',[150 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx2(IndexY,IndexX,1)*len2, Dy2(IndexY,IndexX,1)*len2,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx2(IndexY,IndexX,2)*len2, Dy2(IndexY,IndexX,2)*len2,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Expansion and contraction.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;


% *************************************************************************
% Estimate and display flow for the 3rd sequence.
% *************************************************************************
addpath(algoPath); 
[Dx3 Dy3] = estimateOpticFlow2D(ImgSeq3); 
rmpath(algoPath);
figure('Position',[350 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx3(IndexY,IndexX,1)*len3, Dy3(IndexY,IndexX,1)*len3,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx3(IndexY,IndexX,2)*len3, Dy3(IndexY,IndexX,2)*len3,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Clockwise and counterclockwise rotation.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));