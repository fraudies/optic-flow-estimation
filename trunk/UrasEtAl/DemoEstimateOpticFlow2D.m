clc
clear all
close all

% *************************************************************************
% Uras et al. (1988).
%   Computes and displays the optic flow and the validity of its estimate
%   for expansion.
%   Copyright (C) 2013  Florian Raudies, 01/14/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Generate stimulus.
addpath('./../'); % add root directory to path to access auxiliary functions.
opt.h       = 128; % set parameters for stimulus.
opt.w       = 128;
opt.fNum    = 15;
opt.expRate = 1.005;
[ImgSeq maxSpeed] = randomTextureExpansion(opt);
rmpath('./../');

% Estimate optic flow for the image sequence.
[Dx Dy Valid] = estimateOpticFlow2D(ImgSeq);

% Select first frame.
Dx = squeeze(Dx(:,:,1));
Dy = squeeze(Dy(:,:,1));
Valid = squeeze(Valid(:,:,1));

% Display the estimated optic flow.
h       = size(ImgSeq,1);
w       = size(ImgSeq,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h; 
len     = sample/maxSpeed;
SampleGrid = zeros(h,w);
SampleGrid(IndexY,IndexX) = 1;
% Define a grid of valid sample points.
SampleGrid = SampleGrid > 0;
SampleValid = Valid & SampleGrid;
SampleNotValid = ~Valid & SampleGrid;

figure('Position',[50 50 600 600]); 
quiver(X(SampleValid),      Y(SampleValid),...
       Dx(SampleValid)*len, Dy(SampleValid)*len,0,'-b');
hold on;
quiver(X(SampleNotValid),      Y(SampleNotValid),...
       Dx(SampleNotValid)*len, Dy(SampleNotValid)*len,0,'-r');
legend('valid','invalid');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));

