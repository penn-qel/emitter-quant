global sample_area N blur

%Simulation Size
sample_area = 100; %um^2

%Fixed defects
fixed_density_v = [];
fixed_brightness_v = [];
num_fixed_defect_species = length(fixed_density_v);

%Line defect, goes diagonally through the sample
line_brightness = 0;

%Defect parameters
%Normal parameters
%density_v = [0.03, 0.3, 3, 30];
%brightness_v = [30, 100, 300, 1000];

%Wide Parameters
%density_v = [0.01, 1, 100];
%brightness_v = [10, 100, 1000];

%Single Parameters
density_v = [1];
brightness_v = [100];

%How close initial guesses are
initial_guess_error = 0.1;

%Background parameters
lambda = 100;

%Optical blur
pixel_size = 100; %nm
gaussian_blur = 150; %nm
blur = gaussian_blur / pixel_size;

%Markers for scatter_plot
true_marker = 'o';
fit_marker = 'x';

%Calculating for later
N = round(sqrt(sample_area) / (pixel_size/1000)); %Linear size of system, in pixels

%Seed the rng for reproducibility
rng(0);
figure;
plot_height = 4*length(density_v);
plot_width = 2*length(brightness_v);
pl_map_vector = @(iDensity, iBrightness) [plot_width*(4*iDensity-4)+2*iBrightness-1, plot_width*(4*iDensity-4)+2*iBrightness, plot_width*(4*iDensity-3)+2*iBrightness-1, plot_width*(4*iDensity-3)+2*iBrightness];
histogram_vector = @(iDensity, iBrightness) [plot_width*(4*iDensity-1)+2*iBrightness-1, plot_width*(4*iDensity-1)+2*iBrightness];

%Making placeholders for the figure
figure(1);
figure(2);

for iDensity = 1:length(density_v)
    for iBrightness = 1:length(brightness_v)
        density = density_v(iDensity);
        brightness = brightness_v(iBrightness);
        width = 0.5*brightness;

        %Adding defects
        [intensity, num_defects] = defect_species_intensity(density, brightness, width);
        density = num_defects / (N*pixel_size/1000)^2; %Correcting the true density
        
        %Adding in background
        intensity = intensity + poissrnd(lambda, N, N);
        
        %Adding fixed defects
        for iFixed = 1:num_fixed_defect_species
            fixed_density = fixed_density_v(iFixed);
            fixed_brightness = fixed_brightness_v(iFixed);
            fixed_width = 0.5*fixed_brightness;
            intensity = intensity + defect_species_intensity(fixed_density, fixed_brightness, fixed_width);
        end
        
        %Adding line defect
        intensity = intensity + line_brightness .* imgaussfilt(eye(N), blur);
        
        %% Applying the fitting algorithm
        intensity_vector = intensity(:);
        x0 = [density, brightness, width];
        for iFixed = 1:num_fixed_defect_species
            fixed_density = fixed_density_v(iFixed);
            fixed_brightness = fixed_brightness_v(iFixed);
            fixed_width = 0.5*fixed_brightness;
            x0 = [x0, fixed_density, fixed_brightness, fixed_width];
        end
        x0 = [x0, lambda];
        %x0 = parameters;
        x0 = x0 .* ((1-initial_guess_error) + 2*initial_guess_error*rand(1,4 + 3*num_fixed_defect_species)); %Don't start too close
        W = 1;
        [parameters, err, chi2, hessian, counts, edges] = model_fit_DE(intensity_vector, gaussian_blur, pixel_size, x0, W);
        area = length(intensity_vector);
        AIC = chi2 + 2*length(parameters);
        BIC = chi2 + log(length(counts))*length(parameters);
        [densities, brightnesses, widths, mu, sigma] = extract_params(parameters);
        err_95 = 1.96 * err;
        [density_errs, brightness_errs, width_errs, mu_err, sigma_err] = extract_params(err_95);
        num_species = length(densities);
        
        %% For individual histograms
        figure(1);
        %%Showing the calculated image
        subplot(plot_height, plot_width, pl_map_vector(iDensity ,iBrightness));
        dims = size(intensity);
        xs = pixel_size/1000*(1:dims(1));
        ys = pixel_size/1000*(1:dims(2));
        imagesc(xs, ys, intensity); hold on;
        h = colorbar; set(h, 'YTickLabel', cellstr(num2str(reshape(get(h, 'YTick'),[],1),'%0.0f')) );
        axis square;
        set(gca,'YDir','normal');
        ylabel('(um)');
        xlabel('(um)');
        set(gca,'fontsize', 18);

        
        %% Showing the calculated model
        subplot(plot_height, plot_width, histogram_vector(iDensity ,iBrightness));
        min_x = min(edges) - diff(edges(1:2));
        max_x = max(edges) + diff(edges(1:2));
        histogram('BinEdges',edges,'BinCounts',counts); hold on;
        samples = (edges(1:end-1) + edges(2:end))/2;
        errorbar(samples, counts, sqrt(counts), 'LineStyle', 'none', 'HandleVisibility', 'off');
        x = min_x:sigma/2:max_x;
        plot(x, area.*diff(edges(1:2)).*model_pdf(x, parameters, area, gaussian_blur, pixel_size), 'LineWidth',1, 'Color', 'r');
        %legend('Data', 'Point Defect Model');
        xlabel('Pixel Intensity (Cts)');
        xlim([min(intensity(:)) - diff(edges(1:2)), max(intensity(:)) + diff(edges(1:2))]);
        ylabel('Counts');
        ylim([0.5, Inf]);
        set(gca, 'YScale', 'log');
        %param_string = strcat('FIT PARAMETERS area: ', num2str(area *(pixel_size/1000)^2), ' um^2');
        %param_string = strcat(param_string, '\newline\chi^2: ', num2str(chi2), ' AIC: ', num2str(AIC), ' BIC: ', num2str(BIC));
        param_string = '';
        for i = 1:num_species
            param_string = strcat(param_string, 'Defect Species ', num2str(i), ': \eta= ', num_err(densities(i), density_errs(i)), ' /um^2, I_{avg}= ', num_err(brightnesses(i), brightness_errs(i)), ', \sigma_I= ', num_err(widths(i), width_errs(i)), '\newline');
        end

        if rem(length(x0),3) == 1
            title(strcat(param_string, 'Background: \lambda= ', num_err(mu, mu_err)));

        else
            title(strcat(param_string, 'Background: \mu= ', num_err(mu, mu_err), ', \sigma= ', num_err(sigma, sigma_err)));

        end
        
        
        
        %% For Parameter Plots
        figure(2);
        
        for nSpecies = 1:num_species
            left = brightnesses(nSpecies) - widths(nSpecies);
            right = brightnesses(nSpecies) + widths(nSpecies);
            top = densities(nSpecies)*1.05;
            bottom = densities(nSpecies)/1.05;
            w = fill([left left right right], [bottom, top, top, bottom], 'blue'); hold on;
            w.FaceAlpha = 0.25;
        end
        
        %Showing true densities and brightnesses
        scatter(brightness, density, 100, 'k', true_marker); hold on;
        
        %Showing fitted densities and brightnesses
        scatter(brightnesses, densities, 100, 'b', fit_marker); hold on;
        
        %Fixing errorbars dropping below zero
        density_errs_neg = min(densities - eps(densities), density_errs);
        brightness_errs_neg = min(brightnesses - eps(brightnesses), brightness_errs);
        errorbar(brightnesses, densities, density_errs_neg, density_errs, brightness_errs_neg, brightness_errs, 'b', 'LineStyle', 'none');
        hold on;
        
        %Connect the two
        %plot([brightness brightnesses(1)], [density densities(1)], 'b--'); hold on;
        
        
        %% Showing fixed defect parameters
        figure(2);
        for iFixed=1:num_fixed_defect_species
            fixed_density = fixed_density_v(iFixed);
            fixed_brightness = fixed_brightness_v(iFixed);
            
            %Showing true density and brightness
            scatter(fixed_brightness, fixed_density, 100, true_marker); hold on

            %Connect the two
            plot([fixed_brightness brightnesses(1+iFixed)], [fixed_density densities(1+iFixed)], 'b--'); hold on;
        end
        
        
    end
end
%%
figure(2);
brightness_space = logspace(0, ceil(log10(max(brightness_v))));
background_boundary = 100 * exp(-brightness_space.^2 ./ (2*lambda));
plot(brightness_space, background_boundary, 'r--'); hold on;
plot(xlim(), [(1000/(gaussian_blur))^2 (1000/(gaussian_blur))^2], 'r--'); hold on;
xlim([0.3*min(brightness_v) 3*max(brightness_v)]);
ylim([0.3*min(density_v) 3*max(density_v)]);
xlim([3, 30000]);
ylim([0.03,100]);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('Brightness (cnts)');
ylabel('Density (/um^2)');
title('Emitter Properties');
set(gca,'fontsize', 18);

%% Making a countour plot of the optimization area
model = @(x) area .* diff(edges(1:2)) .* model_pdf(samples, x, area, gaussian_blur, pixel_size);
target = @(x) mod_chi_sq(model(x), counts, 1);

brightness_space = logspace(1, 3, 50);
density_space = logspace(-1, 1, 50);
chi_sq = zeros(length(brightness_space), length(density_space));
for brightness_sample = 1:length(brightness_space)
    for density_sample = 1:length(density_space)
        sample_parameters = parameters;
        sample_parameters(1) = density_space(density_sample);
        sample_parameters(2) = brightness_space(brightness_sample);
        chi_sq(brightness_sample, density_sample) = target(sample_parameters);
    end
end

figure(2);
[Brightness_Samples, Density_Samples] = meshgrid(brightness_space, density_space);
contour(Brightness_Samples, Density_Samples, log10(chi_sq), 10); hold on;
h = colorbar;
ylabel(h, "log(\chi^2)");
%mesh(Brightness_Samples, Density_Samples, log10(chi_sq));
%zlabel("log(\chi^2)")
xlim([10, 1000]);
ylim([0.1,10]);


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

function [intensity, num_defects] = defect_species_intensity(density, brightness, width)
global sample_area N blur

    %Create an array of defects
    num_defects = round(sample_area*density); %Number of defects to place
    brightnesses = brightness + (width/2)*randn(1,num_defects); %Brightnesses
    positions = round(rand(2,num_defects)*N); %Positions
    
    %Create the grid
    x = 0:N-1;
    y = 0:N-1;
    [X,Y] = meshgrid(x,y);
    
    %Placing the defects
    intensity = zeros(N);
    for iDefect = 1:num_defects
        signal = mvnpdf([X(:) Y(:)], positions(:, iDefect)', blur^2*eye(2));
        signal = reshape(signal,length(x),length(y));
        intensity = intensity + brightnesses(iDefect)*signal/max(signal(:));
    end
end