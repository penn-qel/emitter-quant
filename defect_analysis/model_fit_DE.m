function [parameters, err, chi2, hess, counts, edges] = model_fit_DE(intensity, gaussian_blur, pixel_size, x0, W)
	%Estimates model parameters from an intensity map
    
    %Input:
    %intensity = pl map, in counts, unrolled into a vector
    %gaussian_blur = gaussian blurring of the laser, in nm
    %pixel_size = the linear size of pixels, in nm
    
    %x0 = initial parameter guess
    %W = weighting factor
    
    %Output:
    %parameters = model parameters
    %err = estimate of the standard errors for each of the parameters
    %hessian = full hessian of the fit
    %counts = observed counts in each histogram bin
    %edges = edge of histogram bins (nb: length(edges) = length(counts)+1)
    
    % Setting the window
    defect_accuracy = 2; %Window goes from x0./accuracy to x0.*accuracy
    background_accuracy = 1.2;
    
    % How to quantize the brightnesses
    brightness_quantization = 0; %0=Continuous
    
    
    %Deriving some parameters
    area = length(intensity);
    avg = mean(intensity);
    
    [densities, brightnesses, widths, mu, sigma] = extract_params(x0);
    num_species = length(densities);
    
    scale_factor = mu; %Scaling counts by initial background guess
    if rem(length(x0),3) == 1
        scale_vector = [repmat([1, scale_factor, scale_factor],1,num_species), scale_factor];
    else
        scale_vector = [repmat([1, scale_factor, scale_factor],1,num_species), scale_factor, scale_factor];
    end
    x0 = x0 ./ scale_vector;

    %Calculating the counts
    [counts, edges] = histcounts(intensity, ceil(range(intensity)/sqrt(avg)), 'Normalization', 'count');
    
    %Expanding histogram from 0 to twice the maximum intensity
    d = diff(edges(1:2)); %Getting step size
    edges = [fliplr(min(edges):d:0) edges (max(edges):d:2*max(intensity))]; %Being careful to maintain even spacing
    [counts, edges] = histcounts(intensity, edges, 'Normalization', 'count'); %Putting back into histogram function
    
    samples = edges(1:end-1) + diff(edges)/2;

    %Setting up the model target and the starting point
    model = @(x) area .* diff(edges) .* model_pdf(samples, x .* scale_vector, area, gaussian_blur, pixel_size, brightness_quantization);
    target = @(x) mod_chi_sq(model(x), counts, W);
    target_DE = @(num_species, params) mod_chi_sq(model(get_params(num_species, params)), counts, W);
    

    % Set title
    optimInfo.title = 'Histogram Model';
    % Specify objective function
    objFctHandle = target_DE;
    % Define parameter names, ranges and quantization:
    % 1. column: parameter names
    % 2. column: parameter ranges
    % 3. column: parameter quantizations
    % 4. column: initial values (optional)
    paramDefCell = cell(length(x0),4);
    %Defect parameters
    for iSpecies = 1:num_species
        %Density
        paramDefCell{3*iSpecies - 2, 1} = strcat('Density', num2str(iSpecies));
        density_constraint_1 = [0, (1000/pixel_size)^2]; %Between 0 and 1 per pixel
        density_constraint_2 = [(area *(pixel_size/1000)^2)^-1, (1000/(gaussian_blur))^2]; %Between 1 total and 1 per square gaussian blur
        density_accuracy_window = [x0(3*iSpecies - 2) / defect_accuracy, x0(3*iSpecies - 2) * defect_accuracy]; %Assume a window around the initial guess
        paramDefCell{3*iSpecies - 2, 2} = merge_param_ranges(density_constraint_1, density_constraint_2, density_accuracy_window);
        paramDefCell{3*iSpecies - 2, 3} = (area *(pixel_size/1000)^2)^-1 / 10; %Increments of 1 defect
        paramDefCell{3*iSpecies - 2, 4} = x0(3*iSpecies - 2); %Set default to x0
        
        %Brightness
        paramDefCell{3*iSpecies - 1, 1} = strcat('Brightness', num2str(iSpecies));
        brightness_constraint_1 = [0, max(intensity(:))] ./ scale_factor; %Between 0 and the maximum observed brightness
        brightness_constraint_2 = [0, max(intensity(:)) - min(intensity(:))] ./ scale_factor; %Between 0 and the max contrast (assumes one background pixel)
        brightness_accuracy_window = [x0(3*iSpecies - 1) / defect_accuracy, x0(3*iSpecies - 1) * defect_accuracy]; %Assume a window around the initial guess
        paramDefCell{3*iSpecies - 1, 2} = merge_param_ranges(brightness_constraint_1, brightness_constraint_2, brightness_accuracy_window);
        paramDefCell{3*iSpecies - 1, 3} = brightness_quantization ./ scale_factor;
        paramDefCell{3*iSpecies - 1, 4} = x0(3*iSpecies - 1); %Set default to x0
        
        %Width
        paramDefCell{3*iSpecies, 1} = strcat('Width', num2str(iSpecies));
        width_constraint = [0, max(intensity(:))]; %Between 0 and the maximum observed brightness
        width_accuracy_window = [x0(3*iSpecies) / defect_accuracy, x0(3*iSpecies) * defect_accuracy]; %Assume a window around the initial guess
        paramDefCell{3*iSpecies, 2} = merge_param_ranges(width_constraint, width_accuracy_window);
        paramDefCell{3*iSpecies, 3} = brightness_quantization ./ scale_factor;
        paramDefCell{3*iSpecies, 4} = x0(3*iSpecies); %Set default to x0
    end
    
    %Background parameters
    if rem(length(x0),3) == 1
        paramDefCell{end, 1} = 'lambda';
        %lambda_constraint = [min(intensity(:)), max(intensity(:))] ./ scale_factor; %Between the minimum and maximum observed brightness
        lambda_constraint = [0, max(intensity(:))] ./ scale_factor; %Between the minimum and maximum observed brightness
        lambda_accuracy_window = [x0(end) / background_accuracy, x0(end) * background_accuracy]; %Assume a window around the initial guess
        paramDefCell{end, 2} = merge_param_ranges(lambda_constraint, lambda_accuracy_window);
        paramDefCell{end, 3} = brightness_quantization ./ scale_factor; %Increments of 1
        paramDefCell{end, 4} = x0(end); %Set default to x0
    else
        %TODO: Fill in?
        warning('Gaussian Background not yet supported');
    end
    
    %Correcting bounds
    for iParam = 1:length(x0)
        if paramDefCell{iParam, 3}
            paramDefCell{iParam, 2} = round(paramDefCell{iParam, 2} ./ paramDefCell{iParam, 3}) .* paramDefCell{iParam, 3};
        end
    end
    
    % Get default DE parameters
    DEParams = getdefaultparams;
    % Set number of population members (often 10*D is suggested) 
    DEParams.NP = 10*length(x0);
    % Do not use slave processes here. If you want to, set feedSlaveProc to 1 and
    % run startmulticoreslave.m in at least one additional Matlab session.
    DEParams.feedSlaveProc = 0;
    % Set times
    DEParams.maxiter  = 20;
    DEParams.maxtime  = 5*3600; % in seconds
    DEParams.maxclock = [];
    % Set display options
    DEParams.infoIterations = 10;
    DEParams.infoPeriod     = 30; % in seconds
    DEParams.displayResults = false;
    DEParams.saveHistory = false;
    
    %Setting constraint
    DEParams.validChkHandle = @additional_constraints;
    
    % Start differential evolution
    [parameters, chi2, bestFctParams, nrOfIterations, resultFileName] = differentialevolution(...
        DEParams, paramDefCell, objFctHandle, num_species, {}, [], optimInfo);
    parameters = parameters';
    

    %Using the hessian to estimate error
    hess = hessian(@(params) target(params), parameters);
    [eigenvec, eigenval] = eig(hess);
    %TODO: What if there are multiple negative eigenvalues?
    for i = 1:length(parameters)
        if eigenval(i,i) < 0
            warning('Warning - non-positive curvature detected at minimum');
            warning(['Eigenval: ', num2str(eigenval(i,i)), ' Eigenvec: ', num2str(eigenvec(:,i)')]);
            f = @(alpha) target(parameters + alpha*eigenvec(:,i)')-(chi2+1);
            alpha = fzero(f, 0);
            eigenval(i,i) = 2/alpha^2;
        end
    end
    hess = eigenvec*eigenval*inv(eigenvec);
    err = sqrt(diag(inv(hess)));
    err = err';

    err = err .* scale_vector;
    parameters = parameters .*  scale_vector;
    
    function valid = additional_constraints(num_species, params)
        valid = 1;
        x = get_params(num_species, params);
        for iSpecies = 1:num_species
            valid = valid & (((x(3*iSpecies-1) * scale_factor) - (x(3*iSpecies) * scale_factor)^2)<0); %Must be at least poisson width brightness
            valid = valid & ((x(3*iSpecies-1) - x(3*iSpecies))>0); %Brightness > width
            valid = valid & ((x(3*iSpecies-1) + x(end)) * scale_factor < max(intensity(:))); %Brightness + lambda < max intensity
            
            valid = valid &(x(3*iSpecies-2) > (1000/pixel_size)^2 * exp(-(x(3*iSpecies-1) * scale_factor)^2 / (2*scale_factor*x(end)))); %Must rise above noise of poisson background
            
            if (iSpecies > 1)
                %Removing the swap degeneracy in the space TODO: Is this a good idea?
                valid = valid & (x(3*iSpecies-2) >= x(3*(iSpecies-1)-2)); %Ensuring the densities are ordered
            end
        end
    end

    function x = get_params(num_species, params)
        x = zeros(1, 3*num_species + 1);
        for iSpecies = 1:num_species
            x(3*iSpecies - 2) = params.(strcat('Density', num2str(iSpecies)));
            x(3*iSpecies - 1) = params.(strcat('Brightness', num2str(iSpecies)));
            x(3*iSpecies) = params.(strcat('Width', num2str(iSpecies)));
        end
        
        x(end) = params.lambda;
    end

    function param_range = merge_param_ranges(varargin)
        %Merges parameter ranges to the smallest union

        param_range = [-inf inf];
        for iParamRange = 1:nargin
            tmp_freq_range = sort(varargin{iParamRange}); %Make sure they are in ascending order
            param_range(1) = max(param_range(1), tmp_freq_range(1)); %Set the higest minimum
            param_range(2) = min(param_range(2), tmp_freq_range(2)); %Set the lowest maximum
        end
    end
end