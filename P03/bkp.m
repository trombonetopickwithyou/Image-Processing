% Image Segmentation Example

%------------------------------------------------------------------------------------------------
% Demo to illustrate simple blob detection, measurement, and filtering.
% Requires the Image Processing Toolbox (IPT)
% Running time = 7.5 seconds the first run and 2.5 seconds on subsequent runs.
%
% A similar Mathworks demo:
% http://www.mathworks.com/products/image/demos.html?file=/products/demos/shipping/images/ipexprops.html
%
% Code written and posted by ImageAnalyst, July 2009.  Updated April 2015 for MATLAB release R2015a
%------------------------------------------------------------------------------------------------


%% Startup code
close all; clc; clear; clearvars; tic;

fprintf('Starting...\n');
workspace; % Make sure the workspace panel with all the variables is showing.
imtool close all;  % Close all imtool figures.
format long g;
format compact;
captionFontSize = 14;

%% Read in MATLAB demo image (coins)
baseFileName = 'coins.png';
folder = fileparts(which(baseFileName)); % find demo folder
fullFileName = fullfile(folder, baseFileName);

% It doesn't exist in the current folder, look on the search path.
if ~exist(fullFileName, 'file')
	if ~exist(baseFileName, 'file')
		% image not found
		warningMessage = sprintf('Error: the input image file\n%s\nwas not found.\nClick OK to exit the demo.', fullFileName);
		uiwait(warndlg(warningMessage));
		fprintf(1, 'Finished running BlobsDemo.m.\n');
		return;
	end
	% if we found the file, construct the file name.
	fullFileName = baseFileName; % Note: don't prepend the folder.
end

% If we get here, we should have found the image file.
originalImage = imread(fullFileName);

%% Make image grayscale
[rows, columns, numberOfColorChannels] = size(originalImage);

if numberOfColorChannels > 1
	promptMessage = sprintf('Your image file has %d color channels.\nConverting to Grayscale...', numberOfColorChannels);
	originalImage = rgb2gray(originalImage);
end

% display
subplot(2, 2, 1);
imshow(originalImage);
hold on;

set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   % Maximize the figure window.
drawnow;    % Force it to display right now
caption = sprintf('Original image');
title(caption, 'FontSize', captionFontSize);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.


%% Threshold the image
binaryImage = originalImage > 100;          % Bright objects will be chosen (Use '<' for dark objects)
binaryImage = imfill(binaryImage, 'holes'); % get rid of any stray pixels or "holes"

% Display binary image.
subplot(2, 2, 2);
imshow(binaryImage); 
title('image after thresholding', 'FontSize', captionFontSize); 

%% Identify individual blobs

labeledImage = bwlabel(binaryImage, 8);    %makes all pixels in a blob the same value (but unique to other blobs)
blobProps = regionprops(labeledImage, originalImage, 'all'); %extract information calculated during bwlabel()
numberOfBlobs = size(blobProps, 1);

%% get edges
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1);

% plot edges on top of original image
subplot(2, 2, 1);
hold on;
for k = 1:numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end


%% plot blob labels
textFontSize = 14;
labelShiftX = -7;	% align the labels in the centers of the coins.

% center coordinates of ALL the blobs (x and y)
allBlobCentroids = [blobProps.Centroid];
centroidsX = allBlobCentroids(1:2:end-1);
centroidsY = allBlobCentroids(2:2:end);

% Put the labels on the image
subplot(2, 2, 1);
for k = 1 : numberOfBlobs
	text(centroidsX(k) + labelShiftX, centroidsY(k), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold');
end


%% isolate certain blobs using ismember()

allBlobAreas = [blobProps.Area];

allowableAreaIndexes = allBlobAreas < 2000; % Get a list of the blobs that meet our criteria (logical)
keeperIndexes = find(allowableAreaIndexes); % get indices of keeper blobs

% Extract the blobs that meet our criteria (Result will be an image)
keeperBlobsImage = ismember(labeledImage, keeperIndexes);

% Re-label, but only the keeper blobs
labeledDimeImage = bwlabel(keeperBlobsImage, 8);     % Label each blob

% display image with just keepers
subplot(2, 2, 3);
imshow(labeledDimeImage, []);
axis image;
title('Mask of 3 brightest dimes in img', 'FontSize', captionFontSize);


%% mask image to just keeper blobs
maskedImageDime = originalImage;
maskedImageDime(~keeperBlobsImage) = 0;  %set all non-keeper pixels to zero.

%display
subplot(2,2,4);
imshow(maskedImageDime);
axis image;
title('Only the dimes from the original image', 'FontSize', captionFontSize);

%% crop each image individually
elapsedTime = toc;
fprintf('Done making measurements of the features.\n\nElapsed time = %.2f seconds.', elapsedTime);

figure;	% Create a new figure window
set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);   % Maximize the figure window.
    
for k = 1 : numberOfBlobs
    % Find the bounding box of each blob.
    thisBlobsBoundingBox = blobProps(k).BoundingBox;  % Get list of pixels in current blob.

    % Extract each blob into it's own image.
    subImage = imcrop(originalImage, thisBlobsBoundingBox);

    % Display the image.
    subplot(3, 4, k);
    imshow(subImage);
end
