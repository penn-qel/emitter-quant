function cost = mod_chi_sq(model, count, W)
%Optimization cost function - just a modified chi sq
%chi_sq = (model - count).^2 ./ (model + W); %Modified Pearson Chi-sq
chi_sq = (model - count).^2 ./ max(count, W); %Neyman chi-sq
cost = sum(chi_sq);
end

