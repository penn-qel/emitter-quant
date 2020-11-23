function pdf = model_pdf(samples, params, area, gaussian_blur, pixel_size, brightness_quantization)
    %Returns histogram at given samples
    
    %Inputs:
    %samples = data to be returned
    %params = model parameters
    %area = number of pixels
    %gaussian_blur = gaussian width of beam, in nm
    %pixel_size = size of pixels, in nm
    %brightness_quantization = level at which to quantize brightness
    
    %Outputs:
    %pdf = probability density of the intensities given in samples

%Fixing input
if (nargin < 6) || (brightness_quantization == 0)
    brightness_quantization = 1;
end

%Setting some constants
pixel_density = (10^3/pixel_size)^2; %Density of pixels, in 1/um^2    
blur = gaussian_blur / pixel_size; %Gaussian blur, in pixels
radius = sqrt(area); %Size of the area, assuming it's a square

%Extracting parameters
[densities, brightnesses, widths, mu, sigma] = extract_params(params);
num_species = length(densities);
densities = densities ./ pixel_density; %in per pixel units

%Producing an intensity sample space
step = brightness_quantization;
max_I = max(samples(:))+step;
I = 0:step:max_I;
prune = @(a) a(1:length(I));
%Adding in the background
background_pdf = normpdf(I, mu, sigma);
combined_pdf = background_pdf;

%Various species
for i = 1:num_species
    density = densities(i);
    brightness = brightnesses(i);
    width = widths(i);

    %Calculating some other derived parameters
    lambda = area * density; %The expected number of defects in the image

    %Applying the defect brightness model
    pdf1 = single_defect_pdf(I, blur, radius, brightness, width);

    floor_pdf = convpower(pdf1, floor(lambda), step, length(I));
    if floor(lambda) == lambda
        defect_pdf = floor_pdf;
    else
        ceil_pdf = prune(conv(floor_pdf, pdf1)) * step;
        defect_pdf = (ceil(lambda) - lambda)*floor_pdf + (lambda - floor(lambda))*ceil_pdf;
    end
    combined_pdf = prune(conv(combined_pdf, defect_pdf));
end
%Resolving to the output
pdf = interp1(I, combined_pdf, samples(:), 'linear', eps(0));
pdf = reshape(pdf, size(samples));

end