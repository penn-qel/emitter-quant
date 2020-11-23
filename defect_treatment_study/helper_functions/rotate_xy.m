function [X, Y] = rotate_xy(X, Y, theta, center)
%   Registers a cell array of 2d PL scan images
%   adapted from https://www.mathworks.com/matlabcentral/answers/93554-how-can-i-rotate-a-set-of-points-in-a-plane-by-a-certain-angle-about-an-arbitrary-point
    
    %Input:
    %X = x coords of pts
    %Y = y coords of pts
    %theta = rotation angle, in degrees
    %center = center of rotation (assumed to be the same for x and y)
    
    %Output:
    %X = rotated x coords of pts
    %Y = rotated y coords of pts
    theta = theta * pi / 180;
    pts = [X, Y];
    pts = pts.';
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    pts = R*(pts-center) + center;
    pts = pts.';
    X = pts(:,2);
    Y = pts(:,1);
end