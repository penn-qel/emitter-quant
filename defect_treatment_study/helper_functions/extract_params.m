function [densities, brightnesses, widths, mu, sigma] = extract_params(params)
%   Extracts model parameters from feature vector
    
    %Input:
    %parameters = [density1, brightness1, width1, density2, brightness2, width2,
    %..., mu, sigma (if not poisson)]
    
    %Output:
    %densities = [density1, density2, ...]
    %brightnesses = [mu1, mu2, ...]
    %widths = [width1, width2, ...]
    %mu = mu
    %sigma = sigma (if not poisson), sqrt(mu) if poisson
    
if rem(length(params),3) == 1
    mu = params(end);
    sigma = sqrt(mu);
    num_species = (length(params) - 1)/3;
else
    mu = params(end-1);
    sigma = params(end);
    num_species = (length(params) - 2)/3;
end

densities = [];
brightnesses = [];
widths = [];
for i = 1:num_species
    densities = [densities, params(3*i-2)];
    brightnesses = [brightnesses, params(3*i-1)];
    widths = [widths, params(3*i)];
end

end

