function h = my_hist(intensity, limits)
    h = histogram(intensity, 100, 'BinLimits', limits, 'Normalization', 'probability');
    ylabel('Probability');
    xlabel('Photon Counts');
end

