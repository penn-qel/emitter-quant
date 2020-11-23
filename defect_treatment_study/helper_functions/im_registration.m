function [images] = im_registration(images, method, radius)
    %Registers a cell array of 2d PL scan images
    
    %Example:
    %images = im_registration(images);
    
    %images = im_registration(images, 'manual');
    
    %Input:
    %images = cell_array of PL scans with identical parameters
    %method = 'automatic' (uses built in matlab registration) (default)
    % or 'manual'(lets you manually select common reference point)
    %radius = number of pixels for impoint helper circle
    
    %Output:
    %theta = angle of sample relative to straight up and down
if nargin < 2
    method = 'automatic';
end

%% New Style Registration
if strcmp(method,'automatic')
for i = 1:length(images)
    if i ~= 1
        [optimizer, metric] = imregconfig('monomodal');
        tform = imregtform(images{i}, images{i-1}, 'translation', optimizer, metric);
        images{i} = imtranslate(images{i}, [tform.T(3,1), tform.T(3,2)], 'nearest', 'FillValues', 0);
        %images{i} = imwarp(images{i}, tform, 'nearest', 'FillValues', 0);
    end    
end

%% Manual registration
elseif strcmp(method,'manual')
    if nargin < 3
        radius = 40;
    end
    
    callback_fun = @(h) draw_circle(h, radius);
    for i = 1:length(images)
        figure; imagesc(images{i}); grid on; title('Click on center of hole');
        h = impoint;
        addNewPositionCallback(h, callback_fun);
        position = wait(h); close;

        translation = (size(images{i}) ./ 2) - position;
        images{i} = imtranslate(images{i}, translation, 'nearest');
    end

%% Raise error ir not a valid method
else
    error('Not a valid method! Use either `automatic` or `manual`');
end
%{
%% Doing the old style registration
nswps = length(images);
[npts, npts] = size(images{1});
all_pldata = zeros(npts,npts,nswps);
for swp=1:nswps
    all_pldata(:,:,swp) = fliplr(images{swp}');
end

%% Use cross correlation to register all images to the first
refswp = floor(nswps/2); % reference sweep -- all other sweeps will be registered to this one

% interpolate images to higher resolution to resolve sub-pixel offsets
% this helps to eliminate accumulated errors due to rounding
interp_power = 3;
interp_factor = 2^interp_power;
npts_interp = (npts-1)*interp_factor+1;
% images need to be cropped when using cross correlation: this defines the
% number of border pixels to remove 
border = 5*interp_factor; %5
max_drift = 10*interp_factor; % maximum allowed drift (in interpolated pixels) between frames %10*interp_factor

figure;
imshow(images{1}/max(max(images{1})));
message = sprintf('Select a good region for image registration.\nThe algorithm will register based on a bounding rectangle.\nNo selection will select the whole image.');
% Display instructions as title.
title(message, 'FontSize', 16);
uiwait(msgbox(message),60);
% Unfortunately everything disappears when you run roipoly.
[tmp, x, y] = roipoly();
if isempty(tmp) %Default
    crop_rect = [border+1,border+1,npts_interp-2*border-1,npts_interp-2*border-1]; %[xmin ymin width height]
else
    p = polyshape(x,y);
    [xlim, ylim] = boundingbox(p);
    crop_rect = round([xlim(1), ylim(1), xlim(2) - xlim(1), ylim(2) - ylim(1)]*interp_factor); % choose specific sub-region, e.g., to avoid blinking spots  
end

unsort_ixs = 1:nswps;

relative_offsets = zeros(nswps-1,2);
for ss=1:(nswps-1)
    swp1 = ss;
    swp2 = ss+1;
    im1 = interp2(all_pldata(:,:,swp1),interp_power,'linear');
    im2 = interp2(all_pldata(:,:,swp2),interp_power,'linear');
    sub_im2 = imcrop(im2,crop_rect);
    cross_corr = normxcorr2(sub_im2,im1);
    % relative offset of images
    x_center = size(sub_im2,2)+crop_rect(1); % expected index of correlation overlap in x (assuming no drift)
    y_center = size(sub_im2,1)+crop_rect(2); % expected index of " " for y
    allowed_x = x_center+(-max_drift:max_drift);
    allowed_y = y_center+(-max_drift:max_drift);
    corr_mask = zeros(size(cross_corr));
    corr_mask(allowed_y,allowed_x) = 1;
    [~,maxix] = max(cross_corr(:).*corr_mask(:));
    [ypeak,xpeak] = ind2sub(size(cross_corr),maxix);
    offset = [xpeak-x_center, ypeak-y_center]; % this is offset of im1 - im2 
    relative_offsets(ss,:) = -offset/interp_factor; % rescale offset to original pixel size, and minus sign to give offset of im2-im1
end

% Calculate cumulate offsets and reorder to match sweeps
cum_offsets_polorder = [0,0; cumsum(relative_offsets)]; % no offset for 1st image
cum_offsets = cum_offsets_polorder(unsort_ixs,:);

all_offsets = cum_offsets-repmat(cum_offsets(refswp,:),nswps,1); % register cumulative offsets to reference image

% Now register all images to reference sweep
% Now register all images to reference sweep
reg_data_raw = zeros(size(all_pldata));
reg_data_interp = zeros(size(all_pldata));
for ss=1:nswps
    % Register 2nd image to 1st
    swp = ss;
    xform = [1 0 0
        0 1 0
        -all_offsets(ss,1) -all_offsets(ss,2) 1];
    tform = maketform('affine',xform);
    reg_data_raw(:,:,swp) = imtransform(all_pldata(:,:,swp),...
        tform,...
        'nearest',... % No interpolation to preserve raw counts in data for frame-by-frame fitting
        'XData',[1,npts],'YData',[1,npts]);
    reg_data_interp(:,:,swp) = imtransform(all_pldata(:,:,swp),...
        tform,...
        'linear',... % Linear interpolation, better for constructing composite image
        'XData',[1,npts],'YData',[1,npts]);
end


for swp=1:nswps
    images{ss} = reg_data_interp(:,:,ss);
end
%}
end

