%{
%j39
%Info about the data folders
base_data_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\13\j39';
treatment_stages = {'pre_treatment', 'post_sem'};
treatment_stage_strs = {'Pre-Treatment', 'Post-SEM'};
scan = 'Sample 13 j39 polarization';
%Initial parameter Guesses
x0s = {[1.5, 500, 350, 45], [3.3, 650, 600, 100]};
%}

%{
%13, H17
%Info about the data folders
base_data_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\13\H17';
treatment_stages = {'pre_treatment', 'post_sem'};
treatment_stage_strs = {'Pre-Treatment', 'Post-SEM'};
scan = 'Sample 13 H17 polarization';
%Initial parameter Guesses
x0s = {[0.1, 5000, 6000, 1000], [0.1, 5000, 1500, 250]};
%}


%17, k7 8um
title_str = 'Sample 17, Spot k7';
%Info about the data folders
base_data_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\17\k7';
treatment_stages = {'pre_treatment', 'post_anneal', 'post_sem'};
treatment_stage_strs = {'Initial', 'Post-Anneal', 'Post-SEM'};

scan = 'Sample 17 k7 8 um polarization';
%Initial parameter Guesses
x0s = {[0.75, 278, 191, 28.6], [0.78, 1779, 1500, 240], [0.6, 1650, 921, 138]};
%}

%{
%17, k7 4um
title_str = 'Sample 17, Spot k7 4um';
%Info about the data folders
base_data_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\17\k7';
treatment_stages = {'pre_treatment', 'post_anneal', 'post_sem'};
treatment_stage_strs = {'Initial', 'Post-Anneal', 'Post-SEM'};

scan = 'Sample 17 k7 4 um polarization';
%Initial parameter Guesses
x0s = {[1.3, 150, 250, 20], [0.7, 1100, 680, 200], [0.6, 650, 600, 110]};
%}

%{
%A29
title_str = 'Sample 14, Spot A29';
%Info about the data folders
base_data_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\14\A29';
treatment_stages = {'pre_treatment', 'post_anneal', 'post_sem'};
treatment_stage_strs = {'Initial', 'Post-Anneal', 'Post-SEM'};

scan = 'Sample 14 A29 polarization';
%Initial parameter Guesses
x0s = {[0.3, 160, 150, 15], [0.5, 3000, 1000, 100], [0.5, 2000, 800, 100]};
%}

%{
%A35
title_str = 'Sample 17, Spot A35';
%Info about the data folders
base_data_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\17\A35';
treatment_stages = {'pre_treatment', 'post_anneal', 'post_sem'};
treatment_stage_strs = {'Initial', 'Post-Anneal', 'Post-SEM'};

scan = 'Sample 17 A35 polarization';
%Initial parameter Guesses
x0s = {[0.8, 670, 200, 100], [1.0, 3000, 1750, 360], [1.0, 4500, 2000, 100]};
%}

%Number of samples to bootstrap for error
num_bootstrap_samples = 10;

params = [];
errs = [];

for i = 1:length(treatment_stages)
    %Doing the analysis
    data_folder = fullfile(base_data_folder, treatment_stages{i}, scan);
    %delete(fullfile(data_folder, 'mask.mat'));
    [p, e] = defect_analysis_old(data_folder, x0s{i}, num_bootstrap_samples, treatment_stage_strs{i});
    
    %Appending to data
    params = [params; p];
    errs = [errs; e];
end

save(fullfile(base_data_folder, strcat(scan, 'fit parameters')), 'params', 'errs');

%%
figure;
errorbar(1:length(treatment_stages), params(:,1), errs(:,1), '--^','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.75, length(treatment_stages)+0.25]);
xlabel('Treatment stage');
ylabel('Defect Density (1/um^2)');
ylim([0, Inf]);

figure;
errorbar(1:length(treatment_stages), params(:,2), errs(:,2), '--o','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.75, length(treatment_stages)+0.25]);
xlabel('Treatment stage');
ylabel('Defect Brightness (cts)');
ylim([0, Inf]);

figure;
errorbar(1:length(treatment_stages), params(:,3), errs(:,3), '--s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.75, length(treatment_stages)+0.25]);
xlabel('Treatment stage');
ylabel('Background Brightness (cts)');
ylim([0, Inf]);

%{
figure;
errorbar(1:length(treatment_stages), params(:,1),errs(:,1)); hold on; yyaxis right; errorbar(1:length(treatment_stages), params(:,2),errs(:,2));
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.25, length(treatment_stages)+0.25]);

title(strcat('Defect parameters estimates, ', title_str));
xlabel('Treatment stage');
yyaxis left;
ylabel('Defect Density (1/um^2)');
ylim([0, Inf]);
yyaxis right;
ylabel('Defect Brightness (AU)');
ylim([0, Inf]);

figure;
errorbar(1:length(treatment_stages), params(:,3),errs(:,3)); hold on; yyaxis right; errorbar(1:length(treatment_stages), params(:,4),errs(:,4));
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.25, length(treatment_stages)+0.25]);

%title(strcat('Background parameters estimates, ', title_str))
xlabel('Treatment stage');
yyaxis left;
ylabel('Background Mean (AU)');
ylim([0, Inf]);
yyaxis right;
ylabel('Background Std (AU)');
ylim([0, Inf]);
%}