%Constants of our samples
angle1 = -3.2623; %Pre-treatment
angle2 = -0.7463; %Post-sem

%i31
%{
radius = 40;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\i31_8um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\i31_8um_polarization';
%}
%J4
%{
radius = 30;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\J4_6um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\J4_6um_polarization';
%}
%O13
%{
radius = 40;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\O13_8um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\O13_8um_polarization';
%}
%O12
%{
radius = 30;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\O12_6um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\O12_6um_polarization';
%}
%O35

radius = 40;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\O35_8um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\O35_8um_polarization';
%}
%W7
%{
radius = 40;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\W7_8um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\W7_8um_polarization';
%}
%Y17
%{
radius = 40;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\Y17_8um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\Y17_8um_polarization';
%}
%Z27
%{
radius = 20;
data_folder1 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\Z27_4um_polarization';
data_folder2 = 'S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\Z27_4um_polarization';
%}

channel = 1; %Channel index of the counter
sample_pols = 0:15:165;

%Loading the images
[images_1, pols_1] = load_images(data_folder1, channel);
[images_2, pols_2] = load_images(data_folder2, channel);

%%Registering the images
images_1 = im_registration(images_1, 'manual', radius);
images_2 = im_registration(images_2, 'manual', radius);

%%Taking the FFT
im_fft_1 = images_pol_fft(images_1, pols_1, sample_pols);
im_fft_2 = images_pol_fft(images_2, pols_2, sample_pols);

intensity_1 = im_fft_1{1};
intensity_2 = im_fft_2{1};

%%Making same size by zero padding, assuming the images are square!
if size(intensity_1,1) > size(intensity_2,1)
    intensity_2 = padarray(intensity_2, size(intensity_1) - size(intensity_2), 'pre');
elseif size(intensity_2,1) > size(intensity_1,1)
    intensity_1 = padarray(intensity_1, size(intensity_2) - size(intensity_1), 'pre');
end

%%Register
images = {intensity_1, intensity_2};
images = im_registration(images, 'manual', radius);
intensity_1 = images{1};
intensity_2 = images{2};

intensity_1 = imrotate(intensity_1, angle1, 'bilinear', 'crop');
intensity_2 = imrotate(intensity_2, angle2, 'bilinear', 'crop');

figure;
%norm = max([intensity_1(:); intensity_2(:)]);
norm = 200;
imshowpair(intensity_1 ./ norm, intensity_2 ./ norm, 'Scaling', 'none');

hold on;
% get axis limits 
x0 = get(gca,'xlim') ;
y0 = get(gca,'ylim') ;
% draw dummy data to show legend 
scatter(0,0,200,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
hold on
scatter(0,0,200,'s','MarkerEdgeColor','m','MarkerFaceColor','m')
% set the mits 
axis([x0 y0])
axis off %hide axis
legend('Pre treatment', 'Post SEM');

set(gca, 'ydir', 'normal');


%{
figure;
colormap('redblue');
imagesc(intensity_2 - intensity_1);
cl = caxis; mx = max(abs(cl)); caxis([-mx, mx]);
colorbar;
set(gca, 'ydir', 'normal');
%}
