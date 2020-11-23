function [parameters, err] = defect_analysis_old(data_folder, x0, title_str)
%   Estimates model parameters from a folder of old style polarization
%   scans
    
    %Input:
    %data_folder = folder with pl scans
    %x0 = initial_parameter_guess
    
    %Output:
    %params = [density, avg_brightness, mu, sigma]
    %err = estimate of the standard errors for each of the parameters
    

%% Some constant parameters
sample_pols = 0:10:180; %Polarizations to sample in the interpolation
gaussian_blur = 225; %nm
pixel_size = 75; %nm
W = 1;

%%Loading and sorting all the data files
[images, pols] = load_images_old(data_folder);

%%Registering the images
images = im_registration(images);

%%Handling the Rotations
rot_file = fullfile(data_folder, 'rot.mat');
if exist(rot_file, 'file')
    load(rot_file);
else
    sample_rotation = input('Enter the clockwise sample rotation in degrees (0 to not rotate)');
    save(rot_file, 'sample_rotation');
end
for i = 1:length(images)
    images{i} = imrotate(images{i}, sample_rotation, 'nearest', 'crop');   
end

%%Taking the FFT
im_fft = images_pol_fft(images, pols, sample_pols);
save(fullfile(data_folder, 'im_fft'), 'im_fft');

%%Calculating the HSV values
[hsv_img, hue, sat, val] = calc_hsv_from_fft(im_fft);
sat(isnan(sat))=0;
rgb_img = hsv2rgb(hsv_img);
save(fullfile(data_folder, 'rgb'), 'rgb_img');

%%Masking to the Region of Interest
mask_file = fullfile(data_folder, 'mask.mat');
if exist(mask_file, 'file')
    load(mask_file);
else
    figure;
    imagesc(abs(im_fft{1}));
    set(gca,'YDir','normal');
    [mask, maskX, maskY] = roipoly();
    save(mask_file, 'mask', 'maskX', 'maskY');
end

figure;
dims = size(im_fft{1});
imagesc(pixel_size/1000*(1:dims(1)), pixel_size/1000*(1:dims(2)), abs(im_fft{1})); hold on;
h = colorbar; set(h, 'YTickLabel', cellstr(num2str(reshape(get(h, 'YTick'),[],1),'%0.0f')) );
hold on;
fill(pixel_size/1000*maskX, pixel_size/1000*maskY, 'g', 'FaceAlpha', 0, 'EdgeColor', 'r', 'LineWidth', 1);
axis square;
set(gca,'YDir','normal');
%title(strcat(title_str, ' Photoluminescence Intensity Map'));
title(strcat('Photoluminescence Map (', title_str, ')'));
legend('Fit Region', 'Location', 'northwest');
ylabel('(um)');
xlabel('(um)');

%% Applying the fitting algorithm
%%Calculating average brightness
intensity = abs(im_fft{1}) / length(sample_pols) * length(pols);
intensity = intensity .* mask;
intensity = intensity(intensity ~= 0);

caxis([0 max(intensity)]);

%%Peforming the fit
%x0 = parameters;
[parameters, err, hessian, counts, edges] = model_fit(intensity, gaussian_blur, pixel_size, x0, W);
area = length(intensity);
[densities, avg_brightnesses, mu, sigma] = extract_params(parameters);
[density_errs, avg_brightness_errs, mu_err, sigma_err] = extract_params(err);
num_species = length(densities);

%% Showing the calculated model
min_x = mu - sigma*4;
max_x = mu + sigma*4 + sum(avg_brightnesses .* sqrt(area * densities * (pixel_size/1000)^2));
max_x = 1000;
figure('Position', [680   558   560   210]);
histogram('BinEdges',edges,'BinCounts',counts); hold on;
x = min_x:sigma/2:max_x;
plot(x, area.*diff(edges(1:2)).*model_pdf(x, parameters, area, gaussian_blur, pixel_size), 'LineWidth',1);
xlim([min_x max_x]);
legend('Data', 'Point Defect Model');
xlabel('Pixel Intensity (Cts)');
ylabel('Counts');

param_string = strcat(title_str, ' FIT PARAMETERS (', num2str(area *(pixel_size/1000)^2), ' um^2)');
for i = 1:num_species
    param_string = strcat(param_string, '\newline Defect Species ', num2str(i), ': \eta= ', num_err(densities(i), density_errs(i)), ' /um^2, I_{avg}= ', num_err(avg_brightnesses(i), avg_brightness_errs(i)));
end

if rem(length(x0),2) == 1
    title(strcat(param_string, '\newline Background: \lambda= ', num_err(mu, mu_err)));

else
    title(strcat(param_string, '\newline Background: \mu= ', num_err(mu, mu_err), ', \sigma= ', num_err(sigma, sigma_err)));

end

function s = num_err(num, err)
    d = -floor(min([log10(err), log10(num)])) + 1;
    s = strcat(num2str(num), '\pm', num2str(err));
end

end