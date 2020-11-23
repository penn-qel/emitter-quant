function [pl, xCoords, yCoords] = get_pl_scans(dataObject)
% function for quickly retrieving pl scans from a bassett-lab data object
% with an arbitrary number of sweeps. Ui control allows for the stepping
% through of the scans.
%
% Inputs:
%   dataObject - a basset-lab data object with a pl scan structure save din
%   the data field
%
% Outputs:
%   pl - a cell array of 2d pl scans, with each cell corresponding to each
%   sweep
%--------------------------------------------------------------------------

% Pull out some necessary data
nSweeps = dataObject.nSweeps;
scanData = dataObject.data;

% check the validity of the scan data
suppliedFields = fields(scanData);
allowedFields = {'plScan','xCoords', 'yCoords', 'clockRate'};

fieldOverlap = setxor(allowedFields, suppliedFields);
if ~isempty(fieldOverlap)
    error('Improper fields in the data structure for a pl scan, must be ''plScan'', ''xCoords'', ''yCoords'', and ''clockRate'' ')
end

pl = {};
xCoords = {};
yCoords = {};

for i = 1:nSweeps
    pl = [pl scanData(i).plScan];
    xCoords = [xCoords scanData(i).xCoords];
    yCoords = [yCoords scanData(i).yCoords];
end
end
