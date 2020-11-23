%%This code is not very clean due to lots of quick changes, venture forth
%%with caution -Alex

%% All new samples, individual plots
collection_name = 'All Flakes';
base_folder = 'S:\Projects\hBN_defects\Data\S32_071618\';
sample_group_titles = {'30s @ 3kV', %i31 8um
                        '30s @ 3kV', %J4 6um
                        '30s @ 5kV', %O35 8um
                        '5m @ 3kV', %012 6um
                        '5m @ 3kV', %O13 8um
                        'No Exposure', %W7 8um
                        'No Exposure', %Y17 8um
                        'No Exposure'}; %Z27 4um
sample_groups = {{'i31_8um'},
                    {'J4_6um'},
                    {'O35_8um'},
                    {'O12_6um'},
                    {'O13_8um'},
                    {'W7_8um'},
                    {'Y17_8um'},
                    {'Z27_4um'}};
mask_name_groups = {{'mask_sem'}, %i31
                    {'mask_anneal'}, %J4
                    {'mask'}, %O35
                    {'mask_anneal'}, %O12
                    {'mask_sem'}, %O13
                    {'mask'}, %W7
                    {'mask'}, %Y17
                    {'mask'}}; %Z27
stages = {'pre treatment', 'post sem', 'post anneal'};
stage_titles = {'Initial', 'Irradiated', 'Annealed'};
allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

%% New samples that were irradiated
collection_name = 'Irradiated Flakes';
base_folder = 'S:\Projects\hBN_defects\Data\S32_071618\';
sample_group_titles = {'30s @ 3kV', %i31 8um
                        '30s @ 3kV', %J4 6um
                        %'30s @ 5kV', %O35 8um
                        '5m @ 3kV', %012 6um
                        '5m @ 3kV'}; %O13 8um
sample_groups = {{'i31_8um'},
                    {'J4_6um'},
                    %{'O35_8um'},
                    {'O12_6um'},
                    {'O13_8um'}};
mask_name_groups = {{'mask_sem'}, %i31
                    {'mask_anneal'}, %J4
                    %{'mask'}, %O35
                    {'mask_anneal'}, %O12
                    {'mask_sem'}}; %O13
stages = {'pre treatment', 'post sem', 'post anneal'};
stage_titles = {'Initial', 'Irradiated', 'Annealed'};
allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

%% New samples, control flakes
collection_name = 'Control Flakes';
base_folder = 'S:\Projects\hBN_defects\Data\S32_071618\';
sample_group_titles = {'No Exposure', %W7 8um
                        'No Exposure', %Y17 8um
                        'No Exposure'}; %Z27 4um
sample_groups = {{'W7_8um'},
                    {'Y17_8um'},
                    {'Z27_4um'}};
mask_name_groups = {{'mask'}, %W7
                    {'mask'}, %Y17
                    {'mask'}}; %Z27
stages = {'pre treatment', 'post sem', 'post anneal'};
stage_titles = {'Initial', 'Irradiated', 'Annealed'};
allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

%% Old annealed-first flakes
base_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans';
sample_group_titles = {'Region D1', 'Region D2', 'Region D3'};
sample_groups = {{'17_k7_8um'}, {'14_A29'}, {'17_A35'}};
mask_name_groups = {{'mask'}, {'mask'}, {'mask', 'mask_anneal', 'mask_anneal'}};
stages = {'pre treatment', 'post anneal', 'post sem'};
stage_titles = {'Initial', 'Annealed', 'Irradiated'};
allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits


marker_size = 200;
num_families = repmat(0*allowed_params, length(stages), 1);
brightness_quantization = 1;
brightness_samples = 10:brightness_quantization:10000;
brightness_samples = brightness_samples';
brightness_distribution = zeros(length(brightness_samples), length(sample_group_titles), length(stages));
combined_brightness_distribution = zeros(length(brightness_samples), length(stages));
noise_floor = brightness_distribution;
combined_noise_floor = combined_brightness_distribution;
areas = zeros(length(sample_group_titles), length(stages));
total_area = 0*1:length(stages);
density_vectors = cell(length(stages), 1);
brightness_vectors = cell(length(stages), 1);

for iGroup = 1:length(sample_group_titles)
    samples = sample_groups{iGroup};
    mask_names = mask_name_groups{iGroup};
    for nsample = 1:length(samples)
        sample = samples{nsample};
        for nstage = 1:length(stages)
            stage = stages{nstage};
            if length(mask_names) == 1
                mask = mask_names{1};
            else
                mask = mask_names{nstage};
            end
            stage_folder = strcat('PL Scans', {' '}, stage);
            stage_folder = stage_folder{1};
            data_folder = fullfile(base_folder, stage_folder, strcat(sample, '_polarization'));
            final_AIC = inf;
            final_params = [];
            final_errs = [];
            for nParams = allowed_params
                param_file = fullfile(data_folder, strcat(mask,'_params', num2str(nParams), '.mat'));
                if exist(param_file, 'file')
                    load(param_file);
                    if AIC < final_AIC && check_params(parameters, err)
                        final_AIC = AIC; final_params = parameters; final_errs = err;
                    end
                end
            end
            parameters = final_params; err = final_errs; AIC = final_AIC;
            [densities, brightnesses, widths, mu, sigma] = extract_params(parameters);
            [density_errs, brightness_errs, width_errs, mu_err, sigma_err] = extract_params(err);
          
            nSpecies = length(densities);
            num_families(nstage, nSpecies + 1) = num_families(nstage, nSpecies+1) + 1;
            
            area_file = fullfile(data_folder, [mask '_area.mat']);
            if exist(area_file, 'file')
                load(area_file);
            else
                area = input(['Enter the area of fit ' sample ' ' stage ' ' mask ': ']);
                save(area_file, 'area');
            end
            areas(iGroup, nstage) = areas(iGroup, nstage) + area;
            
            for iSpecies = 1:nSpecies
                brightness_distribution(:, iGroup, nstage) = brightness_distribution(:, iGroup, nstage) + densities(iSpecies)*normpdf(brightness_samples, brightnesses(iSpecies), widths(iSpecies)/2);
            end
            noise_floor(:, iGroup, nstage) = 100 * exp(-brightness_samples.^2 ./ (2*mu));
            
            density_vectors{nstage} = [density_vectors{nstage} densities];
            brightness_vectors{nstage} = [brightness_vectors{nstage} brightnesses];
        end
    end
end

%{
for nstage = 1:length(stages)
    figure;
    bar(0:length(allowed_params)-1, num_families(nstage, :));
    title([collection_name ' ' stages{nstage} ' (Average Number of Families ' num2str(avgs(nstage)) ')']);
    xlabel('Number of Emitter Families');
    ylabel('Number of Scans');
end

for nstage = 1:length(stages)
    figure;
    histogram(density_vectors{nstage});
    title([collection_name ' ' stages{nstage} ' Densities']);
    xlabel('Density');
    ylabel('Number of Emitter Families');
end

for nstage = 1:length(stages)
    figure;
    histogram(brightness_vectors{nstage});
    title([collection_name ' ' stages{nstage} ' Brightnesses']);
    xlabel('Brightness');
    ylabel('Number of Emitter Families');
end
%}
total_area = sum(areas,1); %Summing over all samples
colors = jet(length(sample_group_titles));
for nstage = 1:length(stages)
    figure;
    
    %Plotting Defect Distributions
    for iGroup = 1:length(sample_group_titles)
        plot(brightness_samples, brightness_samples .* brightness_distribution(:, iGroup, nstage), 'Color', colors(iGroup, :)); hold on;
        combined_brightness_distribution(:, nstage) = combined_brightness_distribution(:, nstage) + areas(iGroup, nstage)*brightness_distribution(:, iGroup, nstage);
    end
    combined_brightness_distribution(:, nstage) = combined_brightness_distribution(:, nstage) / total_area(nstage);
    plot(brightness_samples, brightness_samples .* combined_brightness_distribution(:, nstage), 'k', 'LineWidth', 2); hold on;
    
    %Ploting noise floor
    yl = ylim;
    for iGroup = 1:length(sample_group_titles)
        plot(brightness_samples, noise_floor(:, iGroup, nstage), '--', 'Color', colors(iGroup, :));
        combined_noise_floor(:, nstage) = combined_noise_floor(:, nstage) + areas(iGroup, nstage)*noise_floor(:, iGroup, nstage);
    end
    combined_noise_floor(:, nstage) = combined_noise_floor(:, nstage) / total_area(nstage);
    %plot(brightness_samples, combined_noise_floor(:, nstage), '--k', 'LineWidth', 2); hold on; %TODO: Does this make sense?
    ylim(yl);
    
    num_defects = sum(combined_brightness_distribution(:,nstage)) * brightness_quantization * total_area(nstage);
    title([collection_name ' ' stages{nstage} ' Brightness distribution (# Emitters: ' num2str(num_defects) ' area: ' num2str(total_area(nstage)) 'um^2)']);
    xlabel('Brightness (cts)');
    ylabel('Brightness * Probability Density (cts * emitters/um^2/cts)');
    set(gca, 'XScale', 'log');
    set(gcf, 'Position', [120   138   1020   420]);
    pbaspect([2 1 1])
    legend([sample_group_titles'; {'Combined Distribution'}]);
    ylim([0 100])
end

figure;
for iGroup = 1:length(sample_group_titles)
        plot(brightness_samples, brightness_samples .* brightness_distribution(:, iGroup, 2)- brightness_samples .* brightness_distribution(:, iGroup, 1), 'Color', colors(iGroup, :)); hold on;
end
plot(brightness_samples, brightness_samples .* combined_brightness_distribution(:,2) - brightness_samples .* combined_brightness_distribution(:,1), 'k', 'LineWidth', 2);
density_change = sum(combined_brightness_distribution(:,2) - combined_brightness_distribution(:,1)) * brightness_quantization;
%Depends on the sample group
%title([collection_name ' Irradiated - Exfoliated Brightness distribution (change in density: ' num2str(density_change) '/um^2)']);
title([collection_name ' Annealed - Exfoliated Brightness distribution (change in density: ' num2str(density_change) '/um^2)']);
xlabel('Brightness');
ylabel('Brightness * Probability Density (cts * emitters/um^2/cts)');
set(gca, 'XScale', 'log');
set(gcf, 'Position', [120   138   1020   420]);
pbaspect([2 1 1])
legend([sample_group_titles'; {'Combined Distribution'}]);
%ylim([-2.5 7.5])
ylim([-50 50]);

figure;
for iGroup = 1:length(sample_group_titles)
        plot(brightness_samples, brightness_samples .* brightness_distribution(:, iGroup, 3)-brightness_samples .* brightness_distribution(:, iGroup, 2), 'Color', colors(iGroup, :)); hold on;
end
plot(brightness_samples, brightness_samples .* combined_brightness_distribution(:,3) - brightness_samples .* combined_brightness_distribution(:,2), 'k', 'LineWidth', 2);
density_change = sum(combined_brightness_distribution(:,3) - combined_brightness_distribution(:,2)) * brightness_quantization;
%Depends on the sample group
%title([collection_name ' Annealed - Irradiated Brightness distribution (change in density: ' num2str(density_change) ' /um^2)']);
title([collection_name ' Irradiated - Annealed Brightness distribution (change in density: ' num2str(density_change) '/um^2)']);
xlabel('Brightness');
ylabel('Brightness * Probability Density (cts * emitters/um^2/cts)');
set(gca, 'XScale', 'log');
set(gcf, 'Position', [120   138   1020   420]);
pbaspect([2 1 1])
legend([sample_group_titles'; {'Combined Distribution'}]);
ylim([-50 50])