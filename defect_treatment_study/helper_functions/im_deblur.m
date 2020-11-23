function [images] = im_deblur(images, hsize, sigma)
%IM_DEBLUR Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    hsize = [10 10];
    sigma = 2;
end

PSF = fspecial('gaussian', hsize, sigma);

for i = 1:length(images)
    images{i} = deconvlucy(images{i}, PSF);
end
    
end

