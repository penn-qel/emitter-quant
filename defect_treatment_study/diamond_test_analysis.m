data_folder = 'S:\Projects\hBN_defects\Papers\TreatmentEffects\Data\diamond_pl_scans';
data_filename = '1bscan002_062719.mat';
data_file = fullfile(data_folder, data_filename);
channel = 2;
gaussian_blur = 210; %nm

%True values
true_brightnesses = [290, 110];
true_densities = [0.0075, 0.025];

%Markers for scatter_plot
true_marker = 'o';
fit_marker = 'x';

load(data_file);
intensity = data.plScan(:, :, channel) * 1000 / data.clockRate;
xs = data.xCoords;
ys = data.yCoords;
pixel_size = diff(xs(1:2))*1000;
pixel_size_um = pixel_size / 1000;
[X, Y] = meshgrid(xs, ys);
blur = gaussian_blur / pixel_size;

%% Displaying PL map
figure;
imagesc(xs, ys, intensity); hold on;
h = colorbar; set(h, 'YTickLabel', cellstr(num2str(reshape(get(h, 'YTick'),[],1),'%0.0f')) );
hold on;
axis square;
set(gca,'YDir','normal');
title('Diamond Photoluminescence Intensity Map');
ylabel('(um)');
xlabel('(um)');


%% Applying the fitting algorithm
intensity_vector = intensity(:);
x0 = [0.0075, 290, 30, 0.025, 120, 20, 0.2, 15, 10, 10.1];
[parameters, err, chi2, hessian, counts, edges] = model_fit_DE(intensity_vector, gaussian_blur, pixel_size, x0, 1);
edges  = edges(2:end);
counts = counts(2:end);
area = length(intensity_vector);
AIC = chi2 + 2*length(parameters);
BIC = chi2 + log(length(counts))*length(parameters);
[densities, brightnesses, widths, mu, sigma] = extract_params(parameters);
[density_errs, brightness_errs, width_errs, mu_err, sigma_err] = extract_params(err);
num_species = length(densities);

%% Showing the calculated model
figure;
min_x = min(edges) - diff(edges(1:2));
max_x = max(edges) + diff(edges(1:2));
histogram('BinEdges',edges,'BinCounts',counts); hold on;
samples = (edges(1:end-1) + edges(2:end))/2;
errorbar(samples, counts, sqrt(counts), 'LineStyle', 'none', 'HandleVisibility', 'off');
x = min_x:sigma/2:max_x;
plot(x, area.*diff(edges(1:2)).*model_pdf(x, parameters, area, gaussian_blur, pixel_size), 'LineWidth',1, 'Color', 'r');
legend('Data', 'Point Defect Model');
xlabel('Pixel Intensity (Cts)');
xlim([min(intensity(:)) - diff(edges(1:2)), max(intensity(:)) + diff(edges(1:2))]);
ylabel('Counts');
ylim([0.5, Inf]);
set(gca, 'YScale', 'log');
param_string = strcat('FIT PARAMETERS area: ', num2str(area *(pixel_size/1000)^2), ' um^2');
param_string = strcat(param_string, '\newline\chi^2: ', num2str(chi2), ' AIC: ', num2str(AIC), ' BIC: ', num2str(BIC), '\newline');
for i = 1:num_species
    param_string = strcat(param_string, 'Defect Species ', num2str(i), ': \eta= ', num_err(densities(i), density_errs(i)), ' /um^2, I_{avg}= ', num_err(brightnesses(i), brightness_errs(i)), ', \sigma_I= ', num_err(widths(i), width_errs(i)), '\newline');
end

if rem(length(x0),3) == 1
    title(strcat(param_string, 'Background: \lambda= ', num_err(mu, mu_err)));
else
    title(strcat(param_string, 'Background: \mu= ', num_err(mu, mu_err), ', \sigma= ', num_err(sigma, sigma_err)));
end

%% For Parameter Plots
[densities, brightnesses, widths, mu, sigma] = extract_params(parameters);
err_95 = 0.89*err;
[density_errs, brightness_errs, width_errs, mu_err, sigma_err] = extract_params(err_95);

figure;
%Showing true densities and brightnesses
scatter(true_brightnesses, true_densities, 100, 'k', true_marker); hold on;
%Showing fitted densities and brightnesses
scatter(brightnesses, densities, 100, 'b', fit_marker); hold on;

%Fixing errorbars dropping below zero
density_errs_neg = min(densities - eps(densities), density_errs);
brightness_errs_neg = min(brightnesses - eps(brightnesses), brightness_errs);
errorbar(brightnesses, densities, density_errs_neg, density_errs, brightness_errs_neg, brightness_errs, 'b', 'LineStyle', 'none');
hold on;

%Connect the two
for iFamily = 1:length(true_brightnesses)
    %plot([true_brightnesses(iFamily) brightnesses(iFamily)], [true_densities(iFamily) densities(iFamily)], 'b--'); hold on;
end

brightness_space = logspace(0, ceil(log10(max(true_brightnesses))));
background_boundary = 100 * exp(-brightness_space.^2 ./ (2*mu));
%Noise floor plotting
plot(brightness_space, background_boundary, 'r--'); hold on;
%Resolution limit plotting
%plot(xlim(), [(1000/(gaussian_blur))^2 (1000/(gaussian_blur))^2], 'r--');
xlim([3, 1000]);
ylim([0.003,1]);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('Brightness (cnts)');
ylabel('Density (/um^2)');
title('Defect Parameters');
set(gca,'fontsize', 18);


function s = num_err(num, err)
    if all([abs(num), abs(err)])
        d = -floor(min([log10(abs(err)), log10(abs(num))]))+1;
        if ~(isnan(d))
            num = round(num, d);
            err = round(err, d);
        end
    end
    s = strcat(num2str(num), '\pm', num2str(err));
end