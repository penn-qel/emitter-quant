function [images, pols, xs, ys] = load_images_old(folder, channel)
%LOAD_IMAGES_OLD loads old style pl images from a folder
%   Assumes the folder has matlab data files, named "XXX.mat", where XXX is
%   the polarization angle, with the variable 'pl' in them
files = dir([folder '\*.mat']);
pols = [];
raw_images = {};
single_photon_count = 0.2;
image_size_um = 10;
xs = {};
ys = {};


for file = files'
    pols = [pols str2num(file.name(1:end-4))];
    load([folder '\' file.name]);
    PL = PL ./ single_photon_count;
    raw_images = [raw_images PL];
    
    % TODO: Correctly extract coordinates?
    xs = [xs linspace(0, image_size_um, size(PL, 1)+1)];
    ys = [ys linspace(0, image_size_um, size(PL, 2)+1)];
end
[pols, s] = sort(pols);
images = raw_images(s);
xs = xs(s);
ys = ys(s);
end

