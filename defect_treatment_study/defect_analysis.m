function [parameters] = defect_analysis(data_folder, channel, x0, title_str, mask_name, registration_method, hole_diameter)
%   Estimates model parameters from a folder of new style polarization
%   scans
    
    %Input:
    %data_folder = folder with pl scans
    %x0 = initial_parameter_guess
    
    %Output:
    %params = [density, avg_brightness, mu, sigma]
    %err = estimate of the standard errors for each of the parameters
    

%% Some constant parameters
gaussian_blur = 150; %nm, based on fitting
W = 1;
mask_filename = strcat(mask_name, '.mat');

%% Loading all the data files and computing related parameters
[images, pols, xs, ys] = load_images_old(data_folder, channel);
tmp = xs{1};
pixel_size = diff(tmp(1:2))*1000;
pixel_size_um = pixel_size / 1000;
dims = size(images{1});
xs = pixel_size_um*(1:dims(1));
ys = pixel_size_um*(1:dims(2));
[X, Y] = meshgrid(xs, ys);
hole_radius_pixels = hole_diameter / pixel_size_um / 2;

%% Registering the images
reg_file = fullfile(data_folder, strcat(registration_method,'_registered.mat'));
if exist(reg_file, 'file')
    load(reg_file);
else
    images = im_registration(images, registration_method, hole_radius_pixels);
    save(reg_file, 'images');
    close;
end

%% Raw counts over all scans
intensity = images{1};
for iImage = 2:length(images)
    intensity = intensity + images{iImage};
end

%% Finding Rotation
idcs = strfind(data_folder,'\');
sample_folder = data_folder(1:idcs(end)-1); 
rot_file = fullfile(sample_folder, 'rot.mat');
if exist(rot_file, 'file')
    load(rot_file);
else
    sample_rotation = input('Enter the clockwise sample rotation in degrees (0 to not rotate)');
    save(rot_file, 'sample_rotation');
end

%% Finding Translation
trans_file = fullfile(data_folder, 'trans.mat');
if exist(trans_file, 'file')
    load(trans_file);
else
    callback_fun = @(h) draw_circle(h, hole_radius_pixels);
    figure; imagesc(intensity); grid on; set(gca,'YDir','normal'); title('Click on center of hole');
    h = impoint;
    addNewPositionCallback(h, callback_fun);
    position = wait(h); close;
    sample_translation = (size(intensity) ./ 2) - position;
    save(trans_file, 'sample_translation');
end

%{
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
%save(fullfile(data_folder, 'im_fft'), 'im_fft');

%% Calculating the HSV values
[hsv_img, hue, sat, val] = calc_hsv_from_fft(im_fft);
rgb_img = hsv2rgb(hsv_img);
save(fullfile(data_folder, 'rgb'), 'rgb_img');
%}

%% Masking to the Region of Interest
idcs = strfind(data_folder,'\');
parent_folder = data_folder(1:idcs(end-1)-1); %e.g. S:\Projects\hBN_defects\Data\S32_071618\
pre_treatment_folder = fullfile(parent_folder, 'PL Scans pre treatment'); %e.g. S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\
flake_folder = fullfile(pre_treatment_folder, data_folder(idcs(end)+1:end)); %e.g. S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\O13_8um_polarization
mask_file = fullfile(flake_folder, mask_filename);
if exist(mask_file, 'file')
    load(mask_file);
    d_size = size(intensity,1) - size(mask,1);
    d_pre = floor(abs(d_size) /2);
    d_post = abs(d_size) - d_pre;
    if d_size > 0
        %Padding to size
        mask = padarray(mask, [d_post d_post], 0, 'post');
        mask = padarray(mask, [d_pre d_pre], 0, 'pre');
    elseif d_size < 0
        %Cutting to size
        mask = mask(d_pre+1:end-d_post, d_pre+1:end-d_post);
    end
else
    figure;
    imagesc(intensity);
    set(gca,'YDir','normal');
    title('SELECT REGION OF INTEREST');
    [mask, maskX, maskY] = roipoly();
    mask = imtranslate(mask, sample_translation); %Centering the mask for later use
    maskX = maskX + sample_translation(1);
    maskY = maskY + sample_translation(2);
    mask = imrotate(mask, sample_rotation, 'nearest', 'crop'); %Rotating to standard orientation
    [maskX, maskY] = rotate_xy(maskX, maskY, sample_rotation, (size(intensity,1) ./ 2));
    save(mask_file, 'mask', 'maskX', 'maskY');
    close;
end
mask = imrotate(mask, -sample_rotation, 'nearest', 'crop');
[maskX, maskY] = rotate_xy(maskX, maskY, -sample_rotation, (size(intensity,1) ./ 2));
mask = imtranslate(mask, -sample_translation);
maskX = maskX - sample_translation(1);
maskY = maskY - sample_translation(2);


%% Buffering the mask
mask_buffer = 3 * gaussian_blur / 1000;

mask_poly = polyshape([pixel_size_um*maskX, pixel_size_um*maskY]);
%Intersecting with the field of view
X_l = [min(X(:)) max(X(:))]; Y_l = [min(Y(:)) max(Y(:))];
fov_poly = polyshape([[X_l(1); X_l(2); X_l(2); X_l(1)], [Y_l(1); Y_l(1); Y_l(2); Y_l(2)]]);
int_poly = intersect(mask_poly, fov_poly);

buffer = polybuffer(int_poly, -mask_buffer);
buffered_region = isinterior(buffer, X(:), Y(:));
buffered_region = reshape(buffered_region, dims);
mask = buffered_region;

%% Applying the mask
intensity_vector = intensity .* mask;
intensity_vector = intensity_vector(intensity_vector ~= 0);

%% Displaying PL map
fig1 = figure(1);
imagesc(xs, ys, intensity); hold on;
h = colorbar; set(h, 'YTickLabel', cellstr(num2str(reshape(get(h, 'YTick'),[],1),'%0.0f')) );
hold on;
fill(pixel_size_um*maskX, pixel_size_um*maskY, 'g', 'FaceAlpha', 0, 'EdgeColor', 'r', 'LineWidth', 1); hold on;
plot(buffer);
axis square;
set(gca,'YDir','normal');
title(strcat(title_str, ' Photoluminescence Intensity Map'));
title(strcat('Photoluminescence Map (', title_str, ')'));
legend('Mask', 'Buffered Region', 'Location', 'northwest');
ylabel('(um)');
xlabel('(um)');
caxis([0 max(intensity_vector)]);
savefig(fullfile(data_folder, strcat(mask_name, '_PL Map')));

%% Making inset figure
%{
fig2 = figure(2);
imagesc(xs, ys, unmasked_intensity);
caxis([0 max(intensity_vector)]);
%ylim([6 7.75]);
%xlim([1.5 3.5]);
ylim([5 7]);
xlim([4 6]);
set(gca,'YDir','normal');
[h_m h_i]=inset(fig1,fig2);
colorbar(h_m);
savefig('PL Map inset');
%}

%% Background Estimation
%{
lambdas = [];
errors = [];
cutoffs = 775:775;
for cutoff = cutoffs
bg_vector = intensity_vector .* (intensity_vector < cutoff) ;
bg_vector = floor(bg_vector(find(bg_vector)));
bg_counts = accumarray(bg_vector, 1);
bg_pdf = bg_counts / length(intensity_vector);
bg_values = 1:(cutoff-1);
figure;
plot(bg_values, bg_pdf, 'o'); hold on;
poiss_model = @(x,xdata) x(1)*poisspdf(xdata, x(2))';
[x,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(poiss_model, [0.5, cutoff], bg_values, bg_pdf);
plot(1:2*x(2), poiss_model(x, 1:2*x(2)));
x
covariance = inv(J.'*J)*var(residual)
lambdas = [lambdas, x];
errors = [errors, sqrt(covariance)];
end
%figure;
%errorbar(cutoffs, lambdas, errors);
%max(lambdas)
%min(lambdas)
%}

%% Applying the fitting algorithm
[parameters, err, chi2, hessian, counts, edges] = model_fit_DE(intensity_vector, gaussian_blur, pixel_size, x0, W);
area = length(intensity_vector);
AIC = chi2 + 2*length(parameters);
BIC = chi2 + log(length(counts))*length(parameters);
[densities, brightnesses, widths, mu, sigma] = extract_params(parameters);
[density_errs, brightness_errs, width_errs, mu_err, sigma_err] = extract_params(err);
num_species = length(densities);


%% Showing the calculated model
min_x = min(intensity_vector) - diff(edges(1:2));
max_x = max(intensity_vector) + diff(edges(1:2));
figure('Position', [680   558   560   210]);
histogram('BinEdges',edges,'BinCounts',counts); hold on;
samples = (edges(1:end-1) + edges(2:end))/2;
errorbar(samples, counts, sqrt(counts), 'LineStyle', 'none', 'HandleVisibility', 'off'); hold on;
x = min_x:sigma/2:max_x;
plot(x, area.*diff(edges(1:2)).*model_pdf(x, parameters, area, gaussian_blur, pixel_size), 'LineWidth',1, 'Color', 'r');
%xlim([min_x max_x]);
legend('Data', 'Point Defect Model');
xlabel('Pixel Intensity (Cts)');
ylabel('Counts');
ylim([0.5, Inf]);
set(gca, 'YScale', 'log');
param_string = strcat(title_str, ' FIT PARAMETERS area: ', num2str(area *(pixel_size_um)^2), ' um^2');
param_string = strcat(param_string, '\newline\chi^2: ', num2str(chi2), ' AIC: ', num2str(AIC), ' BIC: ', num2str(BIC));
for i = 1:num_species
    param_string = strcat(param_string, '\newline Defect Species ', num2str(i), ': \eta= ', num_err(densities(i), density_errs(i)), ' /um^2, I_{avg}= ', num_err(brightnesses(i), brightness_errs(i)), ', \sigma_I= ', num_err(widths(i), width_errs(i)));
end

if rem(length(x0),3) == 1
    title(strcat(param_string, '\newline Background: \lambda= ', num_err(mu, mu_err)));

else
    title(strcat(param_string, '\newline Background: \mu= ', num_err(mu, mu_err), ', \sigma= ', num_err(sigma, sigma_err)));

end

%% Saving if better than previous fits
param_file = fullfile(data_folder, strcat(mask_name, '_params', num2str(length(parameters)), '.mat'));
hist_file = fullfile(data_folder, strcat(mask_name, num2str(length(parameters)), '_Histogram Fit'));
if(exist(param_file, 'file'))
    previous_params = load(param_file);
    if previous_params.AIC >= AIC
        savefig(hist_file);
        save(param_file, 'parameters', 'err', 'chi2', 'AIC', 'BIC', 'intensity');
    else
        temp_param_file = fullfile(data_folder, strcat(mask_name, '_params', num2str(length(parameters)), '_TEMP.mat'));
        temp_hist_file = fullfile(data_folder, strcat(mask_name, num2str(length(parameters)), '_Histogram Fit_TEMP'));
        warning('Fit not better than previous result, previous AIC was:');
        warning(num2str(previous_params.AIC));
        warning('Results saved under temp files');
        warning(temp_param_file);
        warning(temp_hist_file);
        savefig(temp_hist_file);
        save(temp_param_file, 'parameters', 'err', 'chi2', 'AIC', 'BIC', 'intensity');
    end
else
    savefig(hist_file);
    save(param_file, 'parameters', 'err', 'chi2', 'AIC', 'BIC', 'intensity');
end

%% Helper functions
function s = num_err(num, err)
    d = -floor(min([log10(err), log10(num)]))+1;
    if ~(isnan(d))
        num = round(num, d);
        err = round(err, d);
    end
    s = strcat(num2str(num), '\pm', num2str(err));
end

end