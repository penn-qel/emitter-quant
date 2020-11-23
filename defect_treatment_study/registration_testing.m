params = [];
errs = [];

for i = 1:length(sample_pols)
    i
    intensity = images{i} .* mask * length(sample_pols);
    intensity = intensity(intensity ~= 0);
    [p, e] = model_fit(intensity, gaussian_blur, pixel_size, x0, W);
    params = [params; p];
    errs = [errs; e];
end

%%
errorbar(sample_pols, params(:,1),errs(:,1)); hold on; yyaxis right; errorbar(sample_pols, params(:,2),errs(:,2));

title('Parameters estimates for different polarizations');
xlabel('Polarization angle');
yyaxis left;
ylabel('Defect Density (1/um^2)');
yyaxis right;
ylabel('Defect Brightness');