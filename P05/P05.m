%programming notes:
%   need to add some sort of way to avoid duplicates...
%   - maybe go backwards?
%   - use the previous frame's worms and labels, and find the closest worm
%     in the future frame that matches it (then remove it from being picked
%     again?)

close all; clear; clc;


%parameters
filename = "Day 4.mov";
threshval = 220;
minpixsize = 80;

v = VideoReader(filename);

%% read in frame
img_bin = getBinaryImage(v, threshval, minpixsize);

%% get individual objects and their properties
labeledImage = bwlabel(img_bin, 8);                     %uniquely label each blob
blobProps = regionprops(labeledImage, img_bin, 'all');  %extract information about each blob
numBlobs = size(blobProps, 1);

%% label each worm

%algorithm idea:
% - On the first run, use regionprops labeling to label worms.
% - On subsequent runs, use x,y coordinates of last frame's labeled worms to
%   determine the labels of each worm. write down how much each worm's
%   coordinates changed by.
% - If it seems that there are two worms together (maybe if numpixels
%   sudenly increases?), assign two labels to it.

% x,y coordinates of all worm blobs
allCentroids = [blobProps.Centroid];
centroidsX = allCentroids(1:2:end-1);
centroidsY = allCentroids(2:2:end);

worms{1}.xpos = centroidsX;
worms{1}.ypos = centroidsY;

for i=1:numBlobs
    worms{1}.label(i) = i;
end

% plot label on top of each worm
hold on;
imshow(img_bin);
for i=1:numBlobs   
    text(centroidsX(i) + 10, centroidsY(i), num2str(i), 'Color', 'white');
end
title(['Number of Worms: ', num2str(numBlobs)]);

    
index = 2;
while hasFrame(v)
    img_bin = getBinaryImage(v, threshval, minpixsize);     %get frame
    
    labeledImage = bwlabel(img_bin, 8);                     %uniquely label each blob
    blobProps = regionprops(labeledImage, img_bin, 'all');  %extract information about each blob
    numBlobs = size(blobProps, 1);
    
    % x,y coordinates of all worm blobs
    allCentroids = [blobProps.Centroid];
    centroidsX = allCentroids(1:2:end-1);
    centroidsY = allCentroids(2:2:end);
    
    worms{index}.xpos = centroidsX;
    worms{index}.ypos = centroidsY;
    
    % use x,y coordinates to assign the label
    for i=1:numBlobs
        worms{index}.label(i) = findClosestWorm(worms{index}.xpos(i), worms{index}.ypos(i), worms{index - 1}.xpos, worms{index - 1}.ypos, worms{index - 1}.label);
    end
    
    
    % show image with labels
    hold on;
    imshow(img_bin);
    for i=1:numBlobs   
        text(worms{index}.xpos(i) + 10, worms{index}.ypos(i), num2str(worms{index}.label(i)), 'Color', 'white');
    end
    title(['Number of Worms: ', num2str(numBlobs)]);
    pause(0);
    
    index = index + 1;
end







%% functions

function bin = getBinaryImage(v, thresh, minpixval)
    rawFrame = readFrame(v);
    
    % Crop video to show just arena, convert to grayscale
    img = imcrop(rawFrame, [400 100 1600 830]);
    img = rgb2gray(img);
    
    %binarize image
    bin = img > thresh;
    
    %remove connected components smaller than min pix_size
    bin = bwareaopen(bin, minpixval);

end

function label = findClosestWorm(targetX, targetY, oldXs, oldYs, oldLabels)
    
    closest = 1;
    for i=1:numel(oldLabels)
        if dist(targetX, targetY, oldXs(i), oldYs(i)) < dist(targetX, targetY, oldXs(closest), oldYs(closest))
            closest = i;
        end
    end

    label = oldLabels(closest);
    
end

function d = dist(x1, y1, x2, y2)
    
    d = sqrt((x1-x2).^2 + (y1-y2).^2);

end
    
    