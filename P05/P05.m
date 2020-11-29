close all; clear; clc;


%parameters
filename = "Day 4.mov";
threshval = 220;
minpixsize = 80;

v = VideoReader(filename);
while hasFrame(v)
    %% read in frame
    img_bin = getBinaryImage(v, threshval, minpixsize);

    %% get individual objects and their properties
    labeledImage = bwlabel(img_bin, 8);                     %uniquely label each blob
    blobProps = regionprops(labeledImage, img_bin, 'all');  %extract information about each blob
    numBlobs = size(blobProps, 1);

    %% label each worm
    
    %algorithm idea:
    % - On the first run, use regionprops labeling...
    % - On subsequent runs, use x,y coordinates of last frame's labeled worms to
    %   determine the labels of each worm. write down how much each worm's
    %   coordinates changed by.
    % - If it seems that there are two worms together (maybe if numpixels
    %   sudenly increases?), assign two labels to it.
    
    worms{1} = blobProps;

    % x,y coordinates of all worm blobs
    allCentroids = [worms{1}.Centroid];
    centroidsX = allCentroids(1:2:end-1);
    centroidsY = allCentroids(2:2:end);

    % plot label on top of each worm
    hold on;
    imshow(img_bin);
    for i=1:numBlobs   
        text(centroidsX(i) + 10, centroidsY(i), num2str(i), 'Color', 'white');
    end
    title(['Number of Worms: ', num2str(numBlobs)]);
    pause(0);
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
