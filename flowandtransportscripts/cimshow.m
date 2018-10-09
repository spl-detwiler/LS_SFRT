% cimshow is a customized version of imshow
%
% cimshow uses 'jet' colormap
%
% range is option -- will default to [min max] if no range is specified
%
function cimshow(x, range)

switch nargin
    case 2

        imshow(x,range,'colormap',jet)
        impixelinfo
    case 1
        imshow(x,[],'colormap',jet)
        impixelinfo
    case 0
        error('need an input image to display')
end
end

