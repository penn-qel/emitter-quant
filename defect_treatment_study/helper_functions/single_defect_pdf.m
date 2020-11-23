function [pdf] = single_defect_pdf(I, blur, radius, brightness, width)
    %Returns brightness pdf from single defect species
    
    %Inputs:
    %I = intensity samples
    
    %blur = gaussian blurring
    %radius = sample radius
    %brightness = average brightness for distribution
    %width = width of brightness distribution
    
    %Outputs:
    %pdf = probability density of intensity from single defect
    
    N = 5;
    interp_factor = 2;
    sigma = abs(width) / 2;
    step = sigma / interp_factor;
    samples = step:step:N*sigma;
    
    heights = normpdf(samples, 0, sigma);
    for i = length(heights)-1:-1:1
        heights(i) = heights(i) - sum(heights(i+1:end));
    end
    norm = sum(2*heights.*samples);
    heights = heights / norm;
    
    pdf = 0*I;
    for i = 1:length(heights)
        pdf = pdf + heights(i)*(flat_dist_pdf(I, blur, radius, brightness, 2*samples(i)));
    end
    
    norm = sum(pdf) * diff(I(1:2));
    pdf = pdf ./ norm;
end