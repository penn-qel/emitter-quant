%%A sample analysis program

%%Physical parameters
gaussian_blur = 225; %nm

%%Loading the intensity
dataObject = load('1a_R9_8um_062618');
[pl, xCoords, yCoords] = get_pl_scans(dataObject);
x = xCoords{1};
pixel_size = (x(2) - x(1))*1000; %nm
intensity = pl{1};

%Masking to the Region of Interest
%mask = roipoly(intensity ./ max(intensity(:))); 
load('sample_analysis_mask'); %Loading pre-selected mask
intensity = intensity .* mask;
intensity = intensity(intensity ~= 0);

%% Performing the fit
%Model Parameters
x0 = [1.0, 25, 10, 10]; %Initial guess for parameters
W = 10; %Weighting factor

[parameters, err, chi2, hessian, counts, edges] = model_fit(intensity, gaussian_blur, pixel_size, x0, W);
err = abs(err);
area = length(intensity);
[density, avg_brightness, width, mu, sigma] = extract_params(parameters);

%% Showing the calculated image
max_x = mu+avg_brightness*4;
figure(1);
subplot(2,1,1)
h = imagesc(pl{1} .* mask);
axis square;
set(gca,'YDir','normal');
title('Intensity Map');
subplot(2,1,2)
histogram('BinEdges',edges,'BinCounts',counts); hold on;
x = linspace(0, max_x);
plot(x, diff(edges(1:2))*area*model_pdf(x, parameters, area, gaussian_blur, pixel_size));
xlim([0 max_x]);
legend('Data','Model');
title(strcat('Defects: \eta= ', num_err(density, err(1)), ' /um^2, I_{avg}= ', num_err(avg_brightness, err(2)), ', \sigma_I= ', num_err(width, err(3)), ' \newline Background: \lambda= ', num_err(mu, err(4))));

function s = num_err(num, err)
    d = -floor(min([log10(err), log10(num)])) + 1;
    num = round(num, d);
    err = round(err, d);
    s = strcat(num2str(num), '\pm', num2str(err));
end
