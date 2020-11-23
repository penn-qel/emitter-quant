%load('S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\W7_8um_polarization\rgb.mat'); %143nm
%load('S:\Projects\hBN_defects\Data\S32_071618\PL Scans post anneal\i31_8um_polarization\rgb.mat'); %203.5nm
%load('S:\Projects\hBN_defects\Data\S32_071618\PL Scans post anneal\Z27_4um_polarization\rgb.mat'); %183 nm
%load('S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\Y17_8um_polarization\rgb.mat'); %158.5nm, %134nm, %143nm
%load('S:\Projects\hBN_defects\Data\S32_071618\PL Scans pre treatment\O12_6um_polarization\rgb.mat'); %154nm, %179nm, %203nm


%hsv = rgb2hsv(rgb_img);
%unmasked_image = hsv(:,:,3);

load('S:\Projects\hBN_defects\Data\S32_071618\PL Scans post sem\O13_8um_polarization\automatic_registered.mat');
unmasked_image =sum(cat(3,images{:}),3);
figure; imagesc(unmasked_image); set(gca,'YDir','normal'); colorbar;
mask = roipoly();
image = unmasked_image .* mask;
sz = size(image);
x = 1:sz(1);
y = 1:sz(2);
[X, Y] = meshgrid(x,y);
xdata = [X(find(mask)) Y(find(mask))];
intensity_v = image(find(mask));
%Guessing initial parameters
x0 = [max(intensity_v) - min(intensity_v), mean(xdata(:,1)), mean(xdata(:,2)), 1.5, min(intensity_v)];
two_d_gauss = @(p, xdata) p(1)*(2*pi*p(4)^2)*mvnpdf(xdata, [p(2), p(3)], p(4)^2*eye(2)) + p(5);
parameters = lsqcurvefit(two_d_gauss, x0, xdata, intensity_v)
prediction = reshape(two_d_gauss(parameters, [X(:) Y(:)]), length(x), length(y));
figure; imagesc(image); title('PL Scan'); set(gca,'YDir','normal'); colorbar;
figure; imagesc(prediction); title('Fit'); set(gca,'YDir','normal'); colorbar;
diff = abs(image - prediction);
diff = diff / max(diff(:));
figure; imagesc(diff); title('Difference'); set(gca,'YDir','normal'); colorbar;

figure;
plot(prediction(:, round(parameters(2)))); hold on;
plot(image(:, round(parameters(2))), '-x');
legend('Fit', 'Data')