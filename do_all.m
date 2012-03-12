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
I = im2double(imread('Checkerboard.jpg'));

I = checkerboard(8,5,5);
I = imresize(im2double(imread('lena.jpg')), [200,200]);
sigma = 10/256;
I = I + sigma*randn(size(I,1), size(I,2));

sig = 1.5;
I1 = gaussian(I, sig);

lam = 0.01; sig=1;
I2 = perona_malik(I, lam, sig);

I3 = nonlinear_diffusion(I, lam, sig);

sfigure; subplot(221); imshow(I); title('original');
subplot(222); imshow(I1); title('gaussian')
subplot(223); imshow(I2); title('PM')
subplot(224); imshow(I3); title('NL-D')

sfigure; 
subplot(221); imagesc(I); title('original');
subplot(222); imagesc(I1-I); title('gaussian')
subplot(223); imagesc(I2-I); title('PM')
subplot(224); imagesc(I3-I); title('NL-D')
colormap(hot(256))

sigma = 10/256;
h = 12*sigma;
a = 2;
I4 = nlmeans(I, h, a);
