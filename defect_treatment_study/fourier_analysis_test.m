rng('default');
A = 100;
blur = 3;
lambda = 10;

a = 100;
area = a^2;
num_emitters = 20;

x = linspace(0,a);
y = linspace(0,a);

pols = 0:10:180;
images = {};
mu = zeros(length(pols),2);
dipoles = zeros(1,length(pols));
for i = 1:num_emitters
    mu(i,:) = a*rand(1,2);
    dipoles(i) = pi*rand();
end


for pol = pols
[X,Y] = meshgrid(x,y);
sum_I = 0*X(:);
sigma = blur^2 * [1 0; 0 1];
for i = 1:num_emitters
    I = mvnpdf([X(:) Y(:)], mu(i,:), sigma)*cos(pol*pi/180 - dipoles(i))^2;
    sum_I = sum_I + I;
end
sum_I = reshape(sum_I, length(x), length(y));
sum_I = A*2*pi*blur^2*sum_I;
sum_I = sum_I + normrnd(lambda, sqrt(lambda), size(sum_I));
images = [images sum_I];
end

im_fft = images_pol_fft(images, pols, sample_pols);

figure;
subplot(2,2,1);
imagesc(im_fft{1});
colorbar;
subplot(2,2,2);
imagesc(abs(im_fft{2}));
colorbar;
subplot(2,2,3);
imagesc(angle(im_fft{2}));
colorbar;
subplot(2,2,4);
imagesc(abs(im_fft{3}));
colorbar;
