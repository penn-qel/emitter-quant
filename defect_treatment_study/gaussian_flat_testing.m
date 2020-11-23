N = 5;
mu = 0;
sigma = 1;
samples_per_sigma = 2;
step = sigma / samples_per_sigma;
samples = (step:step:N*sigma);
heights = normpdf(samples,0, sigma);
for i = length(heights)-1:-1:1
    heights(i) = heights(i) - sum(heights(i+1:end));
end

x = linspace(mu-5*sigma, mu+5*sigma);
y = linspace(0,0);
for i = 1:length(heights)
    y = y + heights(i)*(abs(x-mu) < step*i);
end
norm = sum(2*heights.*samples);
plot(x, y); hold on;
y = y / norm;
plot(x, y); hold on;
plot(x, normpdf(x, mu, sigma));
legend("Raw Approximation", "Normalized Approx.", "True Normal Dist.");
ylabel("Probability Density");
xlabel("\Delta I = I-A");
title("Approximating Normal Distribution with Uniform Distributions")
