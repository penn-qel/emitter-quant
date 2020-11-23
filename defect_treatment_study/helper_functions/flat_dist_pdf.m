function [pdf] = flat_dist_pdf(I, blur, radius, brightness, width)
    %Returns brightness pdf from flat defect distribution
    
    %Inputs:
    %I = intensity samples
    
    %blur = gaussian blurring
    %radius = sample radius
    %brightness = average brightness for distribution
    %width = width of brightness distribution
    
    %Outputs:
    %pdf = probability density of intensity from single defect
    
    max_brightness = brightness + width/2;
    min_brightness = brightness - width/2;
    min_brightness = max(0, min_brightness);

    pdf_max = zero_dist_pdf(I, blur, radius, max_brightness);
    pdf_min = zero_dist_pdf(I, blur, radius, min_brightness);

    pdf = pdf_max*max_brightness - pdf_min*min_brightness;
    pdf = pdf / (max_brightness - min_brightness);
    
end

function [pdf] = zero_dist_pdf(I, blur, radius, max_brightness)
    %Flat distribution that goes from zero to max_brightness
    pdf = -2*blur^2./(I*radius^4).*(pi*radius^2*(I-max_brightness) + 2*blur^2*(max_brightness+I) + 2*blur^2*max_brightness*log(I/max_brightness));
    pdf(isinf(pdf)) = 0;
    pdf(isnan(pdf)) = 0; %For when max_brightness = 0
    pdf = pdf .* (I > max_brightness*exp(-radius^2/(2*blur^2)));
    pdf = pdf .* (I < max_brightness);
    if max_brightness
        pdf = pdf / max_brightness;
    end
    step = I(2) - I(1);
    pdf(1) = 1/step - sum(pdf);
end

