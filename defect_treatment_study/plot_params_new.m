%%This code is not very clean due to lots of quick changes, venture forth
%%with caution -Alex


%% All new samples, individual plots
base_folder = 'S:\Projects\hBN_defects\Data\S32_071618\';
figure_folder = 'S:\Projects\hBN_defects\Papers\TreatmentEffects\Figures\Data';

sample_group_titles = {'Region A1',
                        'Region A2',
                        'Region E1',
                        'Region B1',
                        'Region B2',
                        'Region C1',
                        'Region C2',
                        'Region C3'};
%                        '30s @ 3kV', %i31 8um
%                         '30s @ 3kV', %J4 6um
%                         '30s @ 5kV', %O35 8um
%                         '5m @ 3kV', %012 6um
%                         '5m @ 3kV', %O13 8um
%                         'No Exposure', %W7 8um
%                         'No Exposure', %Y17 8um
%                         'No Exposure'}; %Z27 4um
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
marker_groups = {{'o'}, %i31
                    {'o'}, %J4
                    {'o'}, %O35
                    {'o'}, %O12
                    {'o'}, %O13
                    {'o'}, %W7
                    {'o'}, %Y17
                    {'o'}}; %Z27
stages = {'pre treatment', 'post sem', 'post anneal'};
stage_titles = {'Initial', 'Irradiated', 'Annealed'};
colors = {'b', 'r', 'g'};
allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

sample_group_titles = {'Region A1',
                        'Region A2',
                        'Region E1',
                        'Region B1',
                        'Region B2',
                        'Region C1',
                        'Region C2',
                        'Region C3'};
%                        '30s @ 3kV', %i31 8um
%                         '30s @ 3kV', %J4 6um
%                         '30s @ 5kV', %O35 8um
%                         '5m @ 3kV', %012 6um
%                         '5m @ 3kV', %O13 8um
%                         'No Exposure', %W7 8um
%                         'No Exposure', %Y17 8um
%                         'No Exposure'}; %Z27 4um
sample_groups = {{'i31_8um'}};
mask_name_groups = {{'mask_sem'}};
marker_groups = {{'o'}};
stages = {'pre treatment', 'post sem', 'post anneal'};
stage_titles = {'Initial', 'Irradiated', 'Annealed'};
colors = {'b', 'r', 'g'};
allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

%% Old 13 j39
% base_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\13';
% sample_group_titles = {'Irradiated 13'};
% sample_groups = {{'j39'}};
% mask_name_groups = {{'mask'}};
% marker_groups = {{'o'}};
% stages = {'pre treatment', 'post sem'};
% stage_titles = {'Initial', 'Irradiated'};
% colors = {'b', 'r'};
% allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

%% Old 14 A29 data
% base_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\14';
% sample_group_titles = {'Region D2'};
% sample_groups = {{'A29'}};
% mask_name_groups = {{'mask', 'mask', 'mask'}};
% marker_groups = {{'o'}};
% stages = {'pre treatment', 'post anneal', 'post sem'};
% stage_titles = {'Initial', 'Annealed', 'Irradiated'};
% colors = {'b', 'g', 'r'};
% allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

%% Old 17 k7 data
% base_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans';
% sample_group_titles = {'Region D1'};
% sample_groups = {{'17_k7_8um'}};
% mask_name_groups = {{'mask', 'mask', 'mask', 'mask_second_anneal'}};
% marker_groups = {{'o'}};
% stages = {'pre treatment', 'post anneal', 'post sem', 'second anneal'};
% stage_titles = {'Initial', 'Annealed', 'Irradiated', 'Annealed Again'};
% colors = {'b', 'g', 'r', 'k'};
% allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits

%% Old 17 A35 data
% base_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans';
% sample_group_titles = {'Region D3'};
% sample_groups = {{'17_A35'}};
% mask_name_groups = {{'mask', 'mask_anneal', 'mask_anneal'}};
% marker_groups = {{'o'}};
% stages = {'pre treatment', 'post anneal', 'post sem'};
% stage_titles = {'Initial', 'Annealed', 'Irradiated'};
% colors = {'b', 'g', 'r'};
% allowed_params = [1, 4, 7, 10, 13]; %Just the poisson fits


marker_size = 200;
for iGroup = 1:length(sample_group_titles)
    samples = sample_groups{iGroup};
    mask_names = mask_name_groups{iGroup};
    markers = marker_groups{iGroup};
    figure;
    legend_datasets = [];
    legend_labels = {};
    %subplot(1,2,iGroup); %I hope we never change the number of groups...
    for nsample = 1:length(samples)
        sample = samples{nsample};
        marker = markers{nsample};
        for nstage = 1:length(stages)
            stage = stages{nstage};
            color = colors{nstage};
            %if length(sample_groups) == 1
            %    mask = mask_names{nstage}; %For more difficult samples
            %else
                mask = mask_names{nsample};
            %end
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
            err_95 = 0.89*err;
            [density_errs, brightness_errs, width_errs, mu_err, sigma_err] = extract_params(err_95);
            
            
            map_fig = openfig(fullfile(data_folder, strcat(mask,'_PL Map')));
            %title([sample_group_titles{iGroup} ' PL Map ' stage_titles{nstage}]);
            title('');
            %Adding scale bar
            plot([1 2], [2 2], 'w-', 'LineWidth', 5); hold on;
            xlabel('');
            ylabel('');
            xticks([]);
            yticks([]);
            legend('Suspended Region', 'Buffer');
            colorbar;
            set(gca,'fontsize', 18);
            region_nsp = strrep(sample_group_titles{iGroup}, ' ', '_');
            %savefig(fullfile(figure_folder, [sample_group_titles{iGroup} ' Map ' stage_titles{nstage}]));
            saveas(map_fig, fullfile(figure_folder, [region_nsp '_Map_' stage_titles{nstage} '.png']));
            close;
            
            
            hist_fig = openfig(fullfile(data_folder, strcat(mask, num2str(length(parameters)), '_Histogram Fit')));
            %title([sample_group_titles{iGroup} ' Histogram and Fit ' stage_titles{nstage}]);
            h = findobj(gcf, 'Type', 'histogram');
            xlim([min(h.BinEdges(find(h.Values))) max(h.BinEdges(find(h.Values)+1))]);
            %savefig(fullfile(figure_folder, [sample_group_titles{iGroup} ' Histogram ' stage_titles{nstage}]));
            title('');
            set(gca,'fontsize', 18);
            saveas(hist_fig, fullfile(figure_folder, [region_nsp '_Histogram_' stage_titles{nstage} '.png']));
            close;
            
            

            %Showing densities and brightnesses
            data = scatter(brightnesses, densities, marker_size, marker, 'MarkerEdgeColor', color); hold on;
            %Fixing errorbars dropping below zero
            density_errs_neg = min(densities - eps(densities), density_errs);
            brightness_errs_neg = min(brightnesses - eps(brightnesses), brightness_errs);
            
            errorbar(brightnesses, densities, density_errs_neg, density_errs, brightness_errs_neg, brightness_errs, 'LineStyle', 'none', 'Color', color);
            hold on;
            
            %Recording for legend
            legend_datasets = [legend_datasets data];
            
%             %Showing widths
%             num_species = length(densities);
%             for i = 1:num_species
%                 x = linspace(brightnesses(i) - widths(i)/2, brightnesses(i) + widths(i)/2);
%                 y = 0*x + densities(i);
%                 line(x, y, 'LineWidth', 2, 'Color', color); hold on;
%             end

            %Showing background line
            if brightnesses
                brightness_space = logspace(0, ceil(log10(max(brightnesses))));
            else
                brightness_space = logspace(0, ceil(log10(mu)));
            end
            background_boundary = 100 * exp(-brightness_space.^2 ./ (2*mu));
            plot(brightness_space, background_boundary, ['--' color]);
            
            if densities
                yl = ylim;
                yl(1) = max(yl(1), 0.01);
                ylim([10^floor(log10(min(yl(1), min(densities)))), 10^ceil(log10(max(yl(2), max(densities))))]);
                xl = xlim;
                xl(1) = max(xl(1), 10);
                xlim([10^floor(log10(min(xl(1), min(brightnesses)))), 10^ceil(log10(max(xl(2), max(brightnesses))))]);
            end
            
            xlim([10, 30000]);
            ylim([0.01, 100]);
        end
    end
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('Brightness (cnts)');
    ylabel('Density (/um^2)');
    %title(sample_group_titles{iGroup}, 'Interpreter', 'none');
    title('');
    set(gca,'fontsize', 18);
    legend(legend_datasets, stage_titles);
    %savefig(fullfile(figure_folder, [sample_group_titles{iGroup} ' parameters']));
    saveas(gcf, fullfile(figure_folder, [region_nsp '_parameters.png']));
end