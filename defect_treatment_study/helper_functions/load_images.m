function [images, pols, xs, ys] = load_images(folder, channel)
%LOAD_IMAGES loads save_data style pl images from a folder
%   Assumes the folder has matlab data files, named "*.mat", with
%   all relevant data fields
%   Returns raw counts, polarizations, and xy coords

 %%Loading and sorting all the data files
    files = dir([folder '\1*.mat']);
    pols = [];
    images = {};
    xs = {};
    ys = {};    
    
    for file = files'
        load([folder '\' file.name]);
        PL = data.plScan(:,:, channel) ./ data.clockRate .* 1000;
        images = [images PL];
        pols = [pols sData(7)];
        xs = [xs data.xCoords];
        ys = [ys data.yCoords];
    end
    [pols, s] = sort(pols);
    images = images(s);
    xs = xs(s);
    ys = ys(s);
end