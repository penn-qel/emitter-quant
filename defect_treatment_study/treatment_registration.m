base_data_folder = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\17\k7';
treatment_stages = {'pre_treatment', 'post_anneal', 'post_sem'};
scan = 'Sample 17 k7 8 um polarization\im_fft.mat';

file1 = fullfile(base_data_folder, treatment_stages{1}, scan);
file2 = fullfile(base_data_folder, treatment_stages{3}, scan);

load(file1); I1 = sqrt(abs(im_fft{1})); I1 = I1/max(I1(:));
load(file2); I2 = sqrt(abs(im_fft{1})); I2 = I2/max(I2(:));
%J = imtranslate(I2, [-5, -10],'FillValues',1);
%imshowpair(I1, I2);

%registrationEstimator(I1, I2)

[movingPoints,fixedPoints] = cpselect(I1, I2, 'Wait', true);
mytform = fitgeotrans(movingPoints, fixedPoints, 'projective');
J = imwarp(I1, mytform);
imshowpair(J,I2,'montage')