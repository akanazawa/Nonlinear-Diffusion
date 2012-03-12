%%%%%%%%%%
% Smooth an image using gaussian filter
%%%%%%%%%%
function [I2] = gaussian(I,sigma)
    h = fspecial('gaussian',ceil(sigma)*6, sigma);
    I2 = imfilter(I, h);
end

