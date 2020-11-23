gaussian_blur = 225;
pixel_size = 75;

A = 100;
lambda = 10;
W = 1;
blur = gaussian_blur/pixel_size;
a = 100;
area = a^2;
num_emitters = 40;
rng('default');

x = linspace(0,a);
y = linspace(0,a);

[X,Y] = meshgrid(x,y);
sum_I = 0*X(:);
sigma = blur^2 * [1 0; 0 1];
for i = 1:num_emitters
    mu = a*rand(1,2);
    I = rand()*mvnpdf([X(:) Y(:)], mu, sigma);
    sum_I = sum_I + I;
end
sum_I = reshape(sum_I, length(x), length(y));
sum_I = A*2*pi*blur^2*sum_I;
sum_I = sum_I + normrnd(lambda, sqrt(lambda), size(sum_I));


x0 = [(num_emitters/a^2)*(1000/pixel_size)^2, A/2, lambda, sqrt(lambda)] .* (rand(1,4))*2;
%x0 = parameters;
[parameters, err, hessian, pdf, edges] = model_fit(sum_I(:), gaussian_blur, pixel_size, x0, W);
density = parameters(1);
avg_brightness = parameters(2);
mu = parameters(3);
sigma = parameters(4);


figure(1);
max_x = mu+avg_brightness*2;
subplot(2,1,1)
h = pcolor(X, Y, sum_I);
colorbar;
axis square;
set(gca,'YDir','normal');
title('Intensity');
subplot(2,1,2)
histogram('BinEdges',edges,'BinCounts',pdf); hold on;
x = linspace(min(edges), max_x);
plot(x, area.*diff(edges(1:2)).*model_pdf(x, density, avg_brightness, mu, sigma, area, gaussian_blur, pixel_size));
xlim([0 max_x]);
legend('Data','Model');
title(strcat('Defects: \eta= ', num_err(density, err(1)), ' /um^2, I_{avg}= ', num_err(avg_brightness, err(2)), ' \newline Background: \mu= ', num_err(mu, err(3)), ', \sigma= ', num_err(sigma, err(4))));

function s = num_err(num, err)
    s = strcat(num2str(num), '\pm', num2str(err));
end