function [] = show_rgb(data_folder, registration_method, pols, background, cutoff)
%   Displays RGB map
    
    %Input:
    %data_folder = folder with pl scans
   

%% Registering the images
reg_file = fullfile(data_folder, strcat(registration_method,'_registered.mat'));
if exist(reg_file, 'file')
    load(reg_file);
else
    %error('Registered images not found');
    hole_radius_pixels = 0;
    images = im_registration(images, registration_method, hole_radius_pixels);
    save(reg_file, 'images');
    close;
end

%% Correcting for background and cutoff
for i = 1:length(pols)
    images{i} = min(images{i}-background/length(pols), 2*cutoff/length(pols));
end

%% Taking the FFT
max_pol = 5*round(max(pols)/5); %Round to the nearest 5 degrees?
pol_step = max_pol / (length(pols)-1);
if max_pol == 180  %Get rid of repeated measurement at 180
    max_pol = max_pol - pol_step;
    pols = pols(1:end-1); 
    images = images(1:end-1);
end
sample_pols = 0:pol_step:max_pol;
im_fft = images_pol_fft(images, pols, sample_pols);

%% Calculating the HSV values
[hsv_img, hue, sat, val] = calc_hsv_from_fft(im_fft);
rgb_img = hsv2rgb(hsv_img);
imshow(rgb_img); set(gca,'YDir','normal');
figure;
imagesc(abs(im_fft{2})); set(gca,'YDir','normal');


end