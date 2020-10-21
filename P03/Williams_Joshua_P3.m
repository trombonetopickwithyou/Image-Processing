% Image Segmentation Example

%------------------------------------------------------------------------------------------------
%
% A similar Mathworks demo:
% http://www.mathworks.com/products/image/demos.html?file=/products/demos/shipping/images/ipexprops.html
%
% Code adapted from: "ImageAnalyst" July 2009.  Updated April 2015 for MATLAB release R2015a
%------------------------------------------------------------------------------------------------


%% Startup code
close all; clc; clear; clearvars; tic;

fprintf('Starting...\n');
tic;
workspace; % Make sure the workspace panel with all the variables is showing.
imtool close all;  % Close all imtool figures.
format long g;
format compact;
captionFontSize = 14;

%% Read in image
path = './cards';
images = dir(fullfile(path,'*.jpg'));   % array of all images
numImages = length(images);
fprintf("Found %d images in %s", numImages, path);

% open image
k = 1;
fprintf("\nOpening image #%d: %s\n", k,images(k).name);
orig_img = imread(fullfile(path,images(k).name));

% resize image if necessary
orig_img = imresize(orig_img, [768 1024]);


%display
subplot(2,3,1);
imshow(orig_img);
title('Original');
set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);   % Maximize the figure window.


%% Make image grayscale
[rows, columns, numberOfColorChannels] = size(orig_img);

% if image isn't grayscale
if numberOfColorChannels > 1
	fprintf('Your image file has %d color channels.\nConverting to Grayscale...', numberOfColorChannels);
	orig_img = rgb2gray(orig_img);
    imshow(orig_img);
end
hold on;


%% straighten image
edges = edgeDetect(orig_img);                           % "clean up" image and get edges
[~, ~, top_coords, ~] = getOuterBbox(edges);            % get outer bounding box of card
%edges = imfill(edges, 'holes'); % get rid of any stray pixels or "holes"

subplot(2,3,2);
imshow(edges);

angle = detOrientation(top_coords);                     % determine orientation of card (align to top edge of the card)
rotated_img = imrotate(orig_img, angle);                % rotate image


%% crop image
edges = edgeDetect(rotated_img);                                                % "clean up" image and get edges
[left_coords, right_coords, top_coords, bottom_coords] = getOuterBbox(edges);   % get outer bounding box of card

bbox_xmin = min(left_coords(:,2));
bbox_xmax = max(right_coords(:,2));
bbox_ymin = min(top_coords(:,1));   %"lowest" y-value, but actually the top of the card
bbox_ymax = max(bottom_coords(:,1));

height = abs(bbox_ymin - bbox_ymax);
width = abs(bbox_xmin - bbox_xmax);

%crop
croppedImage = imcrop(rotated_img, [bbox_xmin, bbox_ymin, width, height]);

%display
subplot(2,3,2);
imshow(croppedImage);
hold on;
title('cropped and aligned img');



%% remove pixels outside of bounding box

% get bounding box of rotated img
edges = edgeDetect(croppedImage);
[left_coords, right_coords, top_coords, bottom_coords] = getOuterBbox(edges);

% threshold image, (binarizes it too)
binaryImage = croppedImage < 130;

[numrows, ~] = size(left_coords);
for i=1:numrows
    binaryImage(left_coords(i,1),1:left_coords(i,2)) = 0;
end

[numrows, ~] = size(right_coords);
for i=1:numrows
    binaryImage(right_coords(i,1),right_coords(i,2):end) = 0;
end

[numrows, ~] = size(top_coords);
for i=1:numrows
    binaryImage(1:top_coords(i,1),top_coords(i,2)) = 0;
end

[numrows, ~] = size(bottom_coords);
for i=1:numrows
    binaryImage(bottom_coords(i,1):end,bottom_coords(i,2)) = 0;
end

% remove stray pixels that weren't removed during intitial threshold stage
binaryImage = bwareaopen(binaryImage,10);   %any group of pixels smaller than 10


% display
subplot(2, 3, 3);
imshow(binaryImage); 
title('binary image');



%% Identify individual blobs
labeledImage = bwlabel(binaryImage, 8);                 %makes all pixels in a blob the same value (but unique to other blobs)
blobProps = regionprops(labeledImage, binaryImage, 'all'); %extract information calculated during bwlabel()
numberOfBlobs = size(blobProps, 1);

fprintf('\nfound %d blobs in image\n',numberOfBlobs);

%% find number and suit blobs

% center coordinates of ALL blobs
allBlobCentroids = [blobProps.Centroid];
centroidsX = allBlobCentroids(1:2:end-1);
centroidsY = allBlobCentroids(2:2:end);

% get blob that's highest up on the screen (y value of card number)
[~, num_Index] = min(centroidsY);

% find blob with closest x value (gives us the card suit)
targetX = centroidsX(num_Index);
tmpXvals = centroidsX;
tmpXvals(num_Index) = NaN;

differences = abs(tmpXvals - targetX);
[~, suit_Index] = min(differences, [], 'omitnan');


% subplot(2, 3, 2);
% text(centroidsX(num_Index), centroidsY(num_Index), '7');
% text(centroidsX(suit_Index), centroidsY(suit_Index), '<3');


% isolate number and suit blobs (result is images)
numberBlob = ismember(labeledImage, num_Index);
suitBlob = ismember(labeledImage, suit_Index);


%% crop each blob (need to do resizing)

% number
num_Bbox = blobProps(num_Index).BoundingBox;  % Get list of pixels in current blob.
num_img = imcrop(numberBlob, num_Bbox);  % Extract each blob into it's own image.
num_img = imresize(num_img, [50 30]);


% suit
suit_Bbox = blobProps(suit_Index).BoundingBox;
suit_img = imcrop(suitBlob, suit_Bbox);
suit_img = imresize(suit_img, [40 30]);




% display images
subplot(2, 3, 4);
imshow(num_img);
title('Card Number');

subplot(2,3,5);
imshow(suit_img);
title('Card Suit');


sec = toc;
fprintf('Program finished running in %.2f seconds\n', sec);


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% functions

%'quantizes' pixel values using threshold. 
%uses convolution to determine 2nd derivative.
%finds zero crossings.
function edges = edgeDetect(img)
    tmp=img;    
    threshold = 180;    %(determined experimentally)
    
    % any pixel above threshold is white, anything below is black
    tmp(tmp(:) > threshold) = 255;
    tmp(tmp(:) ~= 255) = 0;

    % determine second derivative using convolution
    h = [0 1 0; 1 -4 1; 0 1 0];         %"Laplacian Kernel"
    edges = edge(tmp,'zerocross',h);    % convolve and find zerocrossings

    % remove stray pixels that weren't removed during intitial threshold stage
    edges = bwareaopen(edges,25);   %any group of pixels smaller than 10
end



%uses zero-crossing data to determine the coordinates of the 
%outer edges of the playing card.
function [left_coords, right_coords, top_coords, bottom_coords] = getOuterBbox(edges)
    [numrows, numcols] = size(edges);

    left_coords = zeros(numrows,2);
    right_coords = zeros(numrows,2);
    for i=1:numrows
        % "find()" gives every coordinate of a nonzero pixel

        % go by each row and and find "outer" non-zero pixels
        allEdges = find(edges(i,:));

        left = min(allEdges);
        if isempty(left)
            left = NaN;
        end

        right = max(allEdges);
        if isempty(right)
           right = NaN; 
        end

        left_coords(i,:) = [i, left];
        right_coords(i,:) = [i, right];
    end

    %remove NaN elements
    left_coords = left_coords(~isnan( left_coords(:,2) ), :);
    right_coords = right_coords(~isnan( right_coords(:,2) ), :);
    
    
    
    % do the same as above, but now for the top and bottom of the card
    top_coords = zeros(numcols,2);
    bottom_coords = zeros(numcols,2);
    
    for i=1:numcols
        
        %go one column at a time, finding the outermost pixels in each
        allEdges = find(edges(:,i));
        
        top = min(allEdges);
        if(isempty(top))
            top = NaN;
        end
        
        bottom = max(allEdges);
        if(isempty(bottom))
           bottom = NaN ;
        end
        
        top_coords(i,:) = [top i];
        bottom_coords(i,:) = [bottom i];
    end
    
    %remove NaN elements
    top_coords = top_coords( ~isnan(top_coords(:,1)) , :);
    bottom_coords = bottom_coords( ~isnan(bottom_coords(:,1)) , :);


end



% determine the orientation of the card, given the coordinates of the left edge.
function angle = detOrientation(top_coords)
    [numrows, ~] = size(top_coords);

    %determine most top corner
    [row,index] = min(top_coords(:,1)); %min y-value
    most_top = [row top_coords(index,2)];

    %determine other two corners (endpoints of "top_coords" matrix)
    corner1 = top_coords(1,:);
    corner2 = top_coords(numrows,:);
    
    
    %see which corner is closest (telling us which side is the short side)
    if(dist(most_top,corner1) < dist(most_top,corner2))
        %furthest = corner1;
        closest_index = 1;
    else
        %furthest = corner2;
        closest_index = numrows;
    end


    %go to a straight spot in between two corners, determine slope there
    half_distance = round((closest_index - index)./2);
    offset = round(abs(half_distance)/2);

    middle_index = index + half_distance;

    s = slope(top_coords, middle_index, offset);
    angle = atand(s);   %get angle in degrees from slope
    
end


% 2-D distance formula
function d = dist(a,b)  %assumes a and b are matrices where the first value is y (row), second is x (colnum)
    d = sqrt( (a(2)-b(2)).^2 + (a(1)-b(1)).^2 );
end

% determines slope (rise/run) at a specific coordinate (index) between
% 'index - offset' and 'index + offset'.
function s = slope(coords, index, offset)
    point_a = coords(index + offset,:);
    point_b = coords(index - offset,:);

    deltay = point_a(1) - point_b(1);
    deltax = point_a(2) - point_b(2);
    
    %rise over run
    s = deltay./deltax;
end



%% backup code


%get edges

% boundaries = bwboundaries(binaryImage);
% numberOfBoundaries = size(boundaries, 1);
% 
% % plot edges on top of original image
% subplot(1, 3, 2);
% hold on;
% for k = 1:numberOfBoundaries
% 	thisBoundary = boundaries{k};
% 	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
% end