function [im_fft] = images_pol_fft(images, pols, sample_pols)
%IMAGES_FFT Takes the fft of images with respect to their polarizations
%Need evenly spaced measurements from 0 to < 180 degrees

image_array = [];
%Helper functions to unroll and reroll the images
unroll = @(x) x(:);
sz = size(images{1});
reroll = @(x) reshape(x, sz);

for i = 1:length(images)
    image_array = [image_array unroll(images{i})];
end
image_array = interp1(pols, image_array.', sample_pols, 'nearest', 'extrap');

im_fft = {};
image_array_fft = fft(image_array);
for i = 1:(length(images))
    im_fft{i} = reroll(image_array_fft(i,:));
end

end