function rgb_img = pol_images_to_rgb(images, pols)
%POL_IMAGES_TO_RGB Converts aligned images at different polarizations to a
%rgb image
%   Takes images, a cell array, describing the same scan at the
%   polarizations specified in pols
    hue = mod(pols, 180) / 180;
    [m,n] = size(images{1});
    mat = zeros(m,n,length(pols));
    for i = 1:length(images)
        mat(:,:,i) = images{i};
    end
    rgb_img = ApplyColorFilter(mat, hue);

    %Copied, with modifications, from S:\Projects\hBN_defects\Analysis\Defect creation paper figures\j39_SEM.m
    low = nanmin(rgb_img(:));
    high = nanmax(rgb_img(:));
    level_lims = [low high]; % [low,high] mapped to full [0,1] color scale
    gamma = 1; % gamma factor (<1 is brighter, >1 is darker)
    rgb_img = imadjust(rgb_img,level_lims,[0 1],gamma);
end

