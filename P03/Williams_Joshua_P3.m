% ECE5367 - P03: Image Segmentation and Detection Project
% Authors:  Joshua Williams
%           Christian Pennington
%
% Description: 
%
%
% Compatability: MATLAB release R2015a or later
% References: ImageAnalyst, "Image Segmentation Example," MATLAB file exchange
%       (https://www.mathworks.com/matlabcentral/fileexchange/25157-image-segmentation-tutorial)
%------------------------------------------------------------------------------------------------

close all; clear; clc;

%% program parameters
%k = 1;  %image number
path = './images';
imgFormat = '*.tif';



%% find images
fprintf("Loading images...\n");

images = dir(fullfile(path,imgFormat));   % array of all images
numImages = length(images);
fprintf("Found %d images in %s\n", numImages, path);


%% load external matrices (for number & suit recognition later)

%numbers
load('numbers.mat','num1', 'num2', 'num3', 'num4', 'num5', 'num6', 'num7','num8','num9','num10','num11','num12','num13');
numMasks = cat(3, num1, num2, num3, num4, num5, num6, num7, num8, num9, num10, num11, num12, num13);
numNames = ["Ace", "2", "3", "4", "5", "6", "7", "8", "9", "10", "Jack", "Queen", "King"];

%suits
load('suits.mat', 'suitHeart', 'suitDiamond', 'suitClub', 'suitSpade');
suitMasks = cat(3, suitHeart, suitDiamond, suitClub, suitSpade);
suitNames = ["Hearts", "Diamonds", "Clubs", "Spades"];



for k=1:numImages
    tic;
    
    %% open one image
    fprintf("\nOpening image #%d: %s\n", k,images(k).name);
    orig_img = imread(fullfile(path,images(k).name));


    %display
    subplot(2,2,1);
    imshow(orig_img);
    title('original image');
    set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);   % Maximize the figure window.




    %% shrink and make image grayscale
    [rows, cols, numberOfColorChannels] = size(orig_img);

    % if image isn't grayscale
    if numberOfColorChannels > 1
        fprintf('input image has %d color channels, converting to grayscale...\n', numberOfColorChannels);
        orig_img = rgb2gray(orig_img);
    end

    % if image is too big
    if (rows > 768) || (cols > 1024)
        rowsOverflow = rows - 768;
        colsOverflow = cols - 1024;

        if rowsOverflow > colsOverflow
            scale = 768./rows;
        else
            scale = 1024./cols;
        end

        fprintf('input image is %d x %d pixels, shrinking to %d x %d...', rows, cols, scale.*rows, scale.*cols);
        orig_img = imresize(orig_img, scale); 
    end



    %% straighten and mask out card
    justCard = straightenCrop(orig_img);

    %display
    subplot(2,2,2);
    imshow(justCard);
    title('rotated and masked card');



    %% extract number and suit objects from the card
    [allObjects, objectProps] = detectObjs(justCard);
    [num_img, suit_img] = findNumSuit(allObjects, objectProps);


    % display
    subplot(2,2,3);
    imshow(num_img);
    title('Card Number');

    subplot(2,2,4);
    imshow(suit_img);
    title('Card Suit');



    %% determine number and suit

    %make number guess
    index = findClosestMatch(num_img, numMasks);    %returns which element in "nums" is the closest match to "num_img"
    numFinal = numNames(index);

    %make suit guess
    index = findClosestMatch(suit_img, suitMasks);
    suitFinal = suitNames(index);


    %% display final guess
    subplot(2,2,1);
    title(sprintf('Card is a %s of %s', numFinal, suitFinal));
    fprintf('\nCard is a %s of %s\n', numFinal, suitFinal);

    sec = toc;
    fprintf('\nProgram finished running in %.2f seconds...\n', sec);




    %% save (when "training" the program)
    % suitSpade = suit_img;
    % fprintf('Saving suit matrix...\n');
    % save('suits.mat','suitSpade', '-append');
    
    pause;

end

close all;
fprintf('\nNo more images to process\n');

%% functions

function justCard = straightenCrop(orig_img)
    % threshold and cleanup image
    binaryImage = orig_img > 130;          %choose bright objects (Use '<' for dark objects)
    binaryImage = imfill(binaryImage, 'holes'); % fill in holes
    binaryImage = bwareaopen(binaryImage,25); %remove small objects

    % identify separate objects
    labeledImage = bwlabel(binaryImage, 8); %makes all pixels in a blob the same value (but unique to other blobs)
    blobProps = regionprops(labeledImage, binaryImage, 'all'); %extract information calculated during bwlabel()

    % choose largest area object (the card)
    allAreas = [blobProps.Area];
    [~, Cardindex] = max(allAreas);
    cardMask = ismember(labeledImage, Cardindex); %result is image

    % rotate card
    allAngles = [blobProps.Orientation];
    angle = 90 - allAngles(Cardindex);
    rotatedCard = imrotate(orig_img, angle);
    cardMask = imrotate(cardMask, angle);   %also rotate mask to mask out object

    % mask image to just the card
    justCard = rotatedCard;
    justCard(~cardMask) = 255;  %make background the same as card whitespace    
end

function [labeledImage, blobProps] = detectObjs(justCard)
    binaryImage = justCard < 130;               % Bright objects will be chosen (Use '<' for dark objects)
    binaryImage = imfill(binaryImage, 'holes'); % fill in object "holes"
    binaryImage = bwareaopen(binaryImage,5);   %remove any group of pixels smaller than 10


    labeledImage = bwlabel(binaryImage, 8);     %makes all pixels in a blob the same value (but unique to other blobs)
    blobProps = regionprops(labeledImage, binaryImage, 'all'); %extract information calculated during bwlabel()
    
    %numberOfBlobs = size(blobProps, 1);
    %fprintf('\nfound %d blobs on the card\n',numberOfBlobs);
end

function [num_img, suit_img] = findNumSuit(allObjects, objectProps)
    % get coordinates of all blobs
    allCentroids = [objectProps.Centroid];
    centroidsX = allCentroids(1:2:end-1);
    centroidsY = allCentroids(2:2:end);

    % find card number object
    [~, num_Index] = min(centroidsY);  

    % find card suit object
    targetX = centroidsX(num_Index);    %xval of number object
    tmpXvals = centroidsX;
    tmpXvals(num_Index) = NaN;
    differences = abs(tmpXvals - targetX);
    [~, suit_Index] = min(differences, [], 'omitnan');  %object with closest xval


    % isolate number and suit blobs (result is images)
    numberBlob = ismember(allObjects, num_Index);
    suitBlob = ismember(allObjects, suit_Index);


    
    % crop and resize

    % number
    num_Bbox = objectProps(num_Index).BoundingBox;
    num_img = imcrop(numberBlob, num_Bbox);
    num_img = imresize(num_img, [50 30]);


    % suit
    suit_Bbox = objectProps(suit_Index).BoundingBox;
    suit_img = imcrop(suitBlob, suit_Bbox);
    suit_img = imresize(suit_img, [40 30]);
end

function guessIndex = findClosestMatch(img, masks) %returns which element in "nums" is the closest match to "num_img"
    [~, ~, numpages] = size(masks);

    sums = zeros(1, numpages);
    errors = zeros(1, numpages);

    %compare this number to database
    for i=1:numpages
        tmp = img;
        %tmp_inv = ~tmp;
        tmp_mask = masks(:,:,i);

        %mask image, see how many pixels overlap with mask
        tmp(~masks(:,:,i)) = 0;    
        sums(i) = sum(tmp(:));
        
        %see how many black pixels the mask overlaps with
        tmp_mask(tmp) = 0;
        %tmp_inv(~masks(:,:,i)) = 0;
        errors(i) = sum(tmp_mask(:));
        
    end

    %assign a value to represent how well each mask ligned up with the number
    vals = sums - errors;
    [~,guessIndex] = max(vals);
end
