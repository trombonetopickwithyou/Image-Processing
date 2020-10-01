clear; clc;
close all;

%Image Processing P02

%Steps:
%1. determine second derivative of image with respect to x and y
%2. add absolute value of two
%3. find zero crossings
%4. rotate image so this bounding box is upright
%5. crop image to bounding box

path = './images';
%path = input('Enter file path: ', 's');

images = dir(fullfile(path,'*.tif')); % array of all .tiff images in path
numImages = length(images);
fprintf("Found %d images in path\n", numImages);

k=1;
img = imread(fullfile(path,images(k).name));

subplot(1,4,1);
imshow(img);
title('original');

%% "clean up" image by using threshold
threshold = 135;

thresh_img=img;
thresh_img(thresh_img(:) > threshold) = 255;
thresh_img(thresh_img(:) ~= 255) = 0;

subplot(1,4,2);    
imshow(thresh_img);
title('Original image after quantizing');

%% determine second derivative using convolution
h = [0 1 0; 1 -4 1; 0 1 0];     %"Laplacian Kernel"
edges = edge(thresh_img,'zerocross',h);

% display image after zerocross with filtering
subplot(1,4,3);
imagesc(edges);
colormap gray
axis image
title('image after edge detection');


%% make bounding box on boundaries

% bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
% Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.



boundaries = bwboundaries(edges, 8, 'noholes');
numberOfBoundaries = size(boundaries, 1);

subplot(1,4,4);
axis image
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
































%% scratch code %%
% while(threshold < 150)
%     tmp=img;
%     tmp(tmp(:) > threshold) = 255;
%     tmp(tmp(:) ~= 255) = 0;
%     
%     % display image after thresholding
%     subplot(1,2,1);    
%     imshow(tmp);
%     title(num2str(threshold));
%     
%     % determine second derivative using convolution
%     h = [0 1 0; 1 -4 1; 0 1 0];     %"Laplacian Kernel"
%     edges = edge(tmp,'zerocross',h);
%     
%     % display image after zerocross with filtering
%     subplot(1,2,2);
%     imagesc(edges);
%     colormap gray
%     axis image
%     
%    
%     threshold = threshold + 2;
%     pause(1);
% end
