function [parameters, err, chi2, hess, counts, edges] = model_fit(intensity, gaussian_blur, pixel_size, x0, W)
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
    [counts, edges] = histcounts(intensity, 2*ceil(range(intensity)/sqrt(avg)), 'Normalization', 'count');
    
    %Expanding histogram from 0 to twice the maximum intensity
    d = diff(edges(1:2)); %Getting step size
    edges = [fliplr(min(edges):d:0) edges (max(edges):d:2*max(intensity))]; %Being careful to maintain even spacing
    [counts, edges] = histcounts(intensity, edges, 'Normalization', 'count'); %Putting back into histogram function
    
    samples = edges(1:end-1) + diff(edges)/2;

    %Setting up the model target and the starting point
    model = @(x) area .* diff(edges) .* model_pdf(samples, x .* scale_vector, area, gaussian_blur, pixel_size);
    target = @(x) mod_chi_sq(model(x), counts, W);

    %Initial simplex optimization over the parameters
    parameters = x0;
    options = optimset('MaxIter', 0);
    [parameters, chi2] = fminsearch(target, x0, options);
    
    
    % Constrained optimization A*x <= b, c(x) < 0
    A = -diag(ones(length(x0),1));
    b = zeros(length(x0),1); % Everything should be positive
    
    pixel_density = (10^3/pixel_size)^2; %Density of pixels, in 1/um^2
    area_um = area / pixel_density; %Area is um^2
    for nSpecies = 1:num_species
        b(3*nSpecies - 2) = -(1/area_um); % At least one of each defect
        A(3*nSpecies, 3*nSpecies) = 1; A(3*nSpecies, 3*nSpecies-1) = -1; %Width less than brightness
    end
    
    if(rem(length(x0),3) == 2) %Background must be positive definite (>= 1) with finite width
        b(end-1) = -10 ./ scale_factor;
        b(end) = -1 ./ scale_factor;
    else
        b(end) = -10 ./ scale_factor;
    end
    
    options = optimoptions('fmincon','Display', 'iter-detailed', 'MaxIter', 1000, 'Algorithm', 'sqp', 'HonorBounds', true, 'RelLineSrchBnd', 1e-3, 'StepTolerance', 1e-6);
    [parameters, chi2, ~, ~, ~, ~, ~] = fmincon(target, parameters, A, b, [], [], [], [], @sigma_constraint, options);
    % Unconstrained optimization
    %options = optimoptions('fminunc', 'MaxIter', 1, 'MaxFunctionEvaluations', 10);
    %[parameters, chi2, exitflag, output, grad, hessian] = fminunc(target, parameters, options);
    

    %Using the hessian to estimate error
    hess = hessian(target, parameters); %https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation
    [eigenvec, eigenval] = eig(hess);
    %TODO: What if there are multiple negative eigenvalues?
    for i = 1:length(parameters)
        if eigenval(i,i) < 0
            warning('Warning - non-positive curvature detected at minimum');
            warning(['Eigenval: ', num2str(eigenval(i,i)), ' Eigenvec: ', num2str(eigenvec(:,i)')]);
            f = @(alpha) target(parameters + alpha*eigenvec(:,i)')-(chi2+9.70);
            alpha = fzero(f, 0);
            eigenval(i,i) = 2*9.70/alpha^2;
        end
    end
    hess = eigenvec*eigenval*inv(eigenvec);
    err = sqrt(diag(inv(hess/9.70)));
    err = err';

    err = err .* scale_vector;
    parameters = parameters .*  scale_vector;
    
    function [c, ceq] = sigma_constraint(x)
        c = zeros(num_species, 1);
        for iSpecies = 1:num_species
            c(iSpecies) = (x(3*iSpecies-1) * scale_factor) - (x(3*iSpecies) * scale_factor)^2; %Must be at least poisson
        end
        if rem(length(x),3) == 2
            c = [c; 0];
            c(end) = (x(end-1) * scale_factor) - (x(end)*scale_factor)^2; %Must be at least poisson
        end
        ceq = [];
    end
end