%%%%%%%%%%
% Nonlinear Diffusion
% CMSC 828 Image Segmentation Project 1
% 
% Runs 4 different smoothing methods:
%  - gaussian
%  - Perona-Malikn
%  - Non-linear diffusion
%  - Non-local Means
%
% Angjoo Kanazawa March 10th 2012
%%%%%%%%%%
addpath('NL-D')
I = im2double(imread('Checkerboard.jpg'));
I0 = checkerboard(8,5,5);
%I0 = imresize(im2double(imread('lena.jpg')), [200,200]);
sigma = 10/256;
I = I0 + sigma*randn(size(I,1), size(I,2));

sig = 1.5;
I1 = gaussian(I, sig);

lam = 0.01; sig=1;
I2 = perona_malik(I, lam, sig);
I2b = pmc(I, lam, 0.1, 10, 1);

I3 = nonlinear_diffusion(I, lam, sig);
I3b = eed(I, lam, sig, 0.1, 3, 1);

sigma = 10/256;
h = 12*sigma;
a = .2;
I4 = nlmeans(I, h, a);


sfigure; subplot(231); imshow(I0); title('original');
subplot(232); imshow(I); title('with noise');
subplot(233);imshow(I1); title('gaussian')
subplot(234); imshow(I2); title('PM')
subplot(235); imshow(I3); title('NL-D')
subplot(236); imshow(I4); title('NL-means')

sfigure; subplot(231); imagesc(I0); title('original');
subplot(232); imagesc(I); title('with noise');
subplot(233);imagesc(I1-I); title('gaussian')
subplot(234); imagesc(I2-I); title('PM')
subplot(235); imagesc(I3-I); title('NL-D')
subplot(236); imagesc(I4-I); title('NL-means')
colorbar
