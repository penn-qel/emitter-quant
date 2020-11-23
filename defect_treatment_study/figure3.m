%User defined parameters
A35_params_file = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\17\A35\Sample 17 A35 polarizationfit parameters.mat';
k7_params_file = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\17\k7\Sample 17 k7 8 um polarizationfit parameters.mat';
A29_params_file = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\14\A29\Sample 14 A29 polarizationfit parameters.mat';
legend_labels = {'Sample 17, Spot A35', 'Sample 17, Spot k7', 'Sample 14, Spot A29'};
treatment_stages = {'pre_treatment', 'post_anneal', 'post_sem'};
treatment_stage_strs = {'Initial', 'Post-Anneal', 'Post-SEM'};

%Loading all the fit parameters
load(A35_params_file);
A35_params = params;
A35_errs = errs;

load(k7_params_file);
k7_params = params;
k7_errs = errs;

load(A29_params_file);
A29_params = params;
A29_errs = errs;


%Making the figure
figure;
subplot(1,3,1);
errorbar(1:length(treatment_stages), A35_params(:,1), A35_errs(:,1), '--^','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on;
errorbar(1:length(treatment_stages), k7_params(:,1), k7_errs(:,1), '--^','MarkerSize',10,...
    'MarkerEdgeColor','green','MarkerFaceColor','green');
hold on;
errorbar(1:length(treatment_stages), A29_params(:,1), A29_errs(:,1), '--^','MarkerSize',10,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.75, length(treatment_stages)+0.25]);
xlabel('Treatment stage');
ylabel('Defect Density (1/um^2)');
ylim([0, Inf]);
legend(legend_labels, 'Location', 'southeast');




subplot(1,3,2);
errorbar(1:length(treatment_stages), A35_params(:,2), A35_errs(:,2), '--o','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on;
errorbar(1:length(treatment_stages), k7_params(:,2), k7_errs(:,2), '--o','MarkerSize',10,...
    'MarkerEdgeColor','green','MarkerFaceColor','green');
hold on;
errorbar(1:length(treatment_stages), A29_params(:,2), A29_errs(:,2), '--o','MarkerSize',10,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.75, length(treatment_stages)+0.25]);
xlabel('Treatment stage');
ylabel('Defect Brightness (counts)');
ylim([0, Inf]);
legend(legend_labels, 'Location', 'southeast');




subplot(1,3,3);
errorbar(1:length(treatment_stages), A35_params(:,3), A35_errs(:,3), '--s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on;
errorbar(1:length(treatment_stages), k7_params(:,3), k7_errs(:,3), '--s','MarkerSize',10,...
    'MarkerEdgeColor','green','MarkerFaceColor','green');
hold on;
errorbar(1:length(treatment_stages), A29_params(:,3), A29_errs(:,3), '--s','MarkerSize',10,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
xticks(1:length(treatment_stages));
xticklabels(treatment_stage_strs);
xlim([0.75, length(treatment_stages)+0.25]);
xlabel('Treatment stage');
ylabel('Background Brightness (counts)');
ylim([0, Inf]);
legend(legend_labels, 'Location', 'southeast');