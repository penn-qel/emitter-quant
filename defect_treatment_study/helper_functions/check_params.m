function valid = check_params(params, errs)
%   Checks if parameters are valid according to study constraints
    
    %Input:
    %params = [density1, brightness1, width1, density2, brightness2, width2,
    %..., mu, sigma (if not poisson)]
    %errs = elementwise error for each parameter
    
    %Output:
    %valid = validity of parameter vector
    
    valid = 1;
    
    [densities, brightnesses, widths, mu, sigma] = extract_params(params);
    [density_errs, brightness_errs, width_errs, mu_err, sigma_err] = extract_params(errs);

    
    if any(density_errs > 2*densities) || any(brightness_errs > 2*brightnesses)
        valid = 0;
    end
    
end

