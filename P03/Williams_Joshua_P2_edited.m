%Joshua Williams
% ECE 5367: Image Processing
% Project 2:  Automatically localize and undo the rotation of a playing card that is presented to a camera (.tif input).

close all; clear; clc;

path = './images';
%path = input('Enter file path (program looks in specified directory for all .tif images): ', 's');



%% load images
images = dir(fullfile(path,'*.tif'));   % array of all .tif images
numImages = length(images);
fprintf("Found %d images in path\n\n", numImages);


%% main code
for k=1:numImages
    clf     %clear figure between images
    
    % open and display unaltered image
    fprintf("Opening image #%d: %s",k,images(k).name);
    img = imread(fullfile(path,images(k).name));
    subplot(1,3,1);
    imshow(img);
    title('Original');
    hold on;
    

    % "clean up" image and get edges
    edges = edgeDetect(img);
    
    % get outer bounding box of card
    [~, ~, top_coords, ~] = getOuterBbox(edges);


    % determine orientation of card (align to top edge of the card)
    angle = detOrientation(top_coords);
    
    % rotate image
    rotated_img = imrotate(img, angle);     %rotate image

    % get bounding box from rotated edge data
    edges = edgeDetect(rotated_img);
    [left_coords, right_coords, top_coords, bottom_coords] = getOuterBbox(edges);
    
    bbox_xmin = min(left_coords(:,2));
    bbox_xmax = max(right_coords(:,2));
    bbox_ymin = min(top_coords(:,1));   %"lowest" y-value, but actually the top of the card
    bbox_ymax = max(bottom_coords(:,1));

    height = abs(bbox_ymin - bbox_ymax);
    width = abs(bbox_xmin - bbox_xmax);

    % crop image
    croppedImage = imcrop(rotated_img, [bbox_xmin, bbox_ymin, width, height]);
    
    
    
    
    %% determine suit and rank
    
    % array of coordinates of the entire edge of the card
    allCoords = [left_coords; right_coords; top_coords; bottom_coords];
    [numrows, numcols] = size(allCoords);
    
    % remove bounding box
    for i=1:numrows
        edges(allCoords(i,1),allCoords(i,2)) = 0;
    end
    edges = bwareaopen(edges,10);   %remove stray groups of pixels
    
    
    % find leftmost pixel
    [~, ~, top_coords, ~] = getOuterBbox(edges);
    
    %determine most topleft corner
    [row,index] = min(top_coords(:,1)); %min col-value
    left = [row top_coords(index,2)];
    
    
    subplot(1,3,3);
    hold on;
    imshow(edges);
    scatter(left(2), left(1));
    
    
    
    
    
    subplot(1,3,2);
    imshow(croppedImage);
    hold on;
    title('Output - Cropped and Aligned');

    fprintf("\nFinished with image. Press any key to continue...\n\n");
    pause;
    
end %for k=1:numImages

disp("Program ended");









%% functions

%'quantizes' pixel values using threshold. 
%uses convolution to determine 2nd derivative.
%finds zero crossings.
function edges = edgeDetect(img)
    tmp=img;    
    threshold = 130;    %(determined experimentally)
    
    % any pixel above threshold is white, anything below is black
    tmp(tmp(:) > threshold) = 255;
    tmp(tmp(:) ~= 255) = 0;

    % determine second derivative using convolution
    h = [0 1 0; 1 -4 1; 0 1 0];         %"Laplacian Kernel"
    edges = edge(tmp,'zerocross',h);    % convolve and find zerocrossings

    % remove stray pixels that weren't removed during intitial threshold stage
    edges = bwareaopen(edges,10);   %any group of pixels smaller than 10
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