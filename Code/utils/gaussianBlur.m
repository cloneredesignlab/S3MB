function result = gaussianBlur(img, filterSize)
    kernel = fspecial('gaussian', filterSize, floor((filterSize+3)/4));
    rows = size(img, 1);
    cols = size(img, 2);
    filterSizeHalf = floor(filterSize/2);
    imgwork = padarray(img, [filterSizeHalf filterSizeHalf], 'replicate');
    result = zeros(rows, cols, 3);
    for i = filterSizeHalf+1:rows+filterSizeHalf
        for j = filterSizeHalf+1:cols+filterSizeHalf
            filter = zeros(1,3);
            for fx = 1:filterSize
                for fy = 1:filterSize
                    filter(1) = filter(1) + kernel(fx, fy) * double(imgwork(i - filterSizeHalf + fx - 1, j - filterSizeHalf + fy - 1, 1));
                    filter(2) = filter(2) + kernel(fx, fy) * double(imgwork(i - filterSizeHalf + fx - 1, j - filterSizeHalf + fy - 1, 2));
                    filter(3) = filter(3) + kernel(fx, fy) * double(imgwork(i - filterSizeHalf + fx - 1, j - filterSizeHalf + fy - 1, 3));
                end
            end
            result(i-filterSizeHalf, j-filterSizeHalf, 1) = filter(1);
            result(i-filterSizeHalf, j-filterSizeHalf, 2) = filter(2);
            result(i-filterSizeHalf, j-filterSizeHalf, 3) = filter(3);
        end
    end
    result = uint8(result);
end




function newimg = padarray(img,padsize)
% padsize = floor((size(img, 2) - size(img, 1)) / 2);
pad = zeros(padsize(1), padsize(2));
extraLineWhenNeeded = zeros(mod(size(img, 2) - size(img, 1), 2) == 1 ,size(img,2)); % Note that extra line will have 0 rows if it's not needed i.e. when the difference between the number of rows and columns of img is even
newimg = [pad; img; pad; extraLineWhenNeeded];
end

% 2D Gaussian filter
function h = gaussian2D(siz, std)

% create the grid of (x,y) values
siz = (siz-1)./2;
[x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));

% analytic function
h = exp(-(x.*x + y.*y)/(2*std*std));

% truncate very small values to zero
h(h<eps*max(h(:))) = 0;

% normalize filter to unit L1 energy
sumh = sum(h(:));
if sumh ~= 0
    h = h/sumh;
end

end