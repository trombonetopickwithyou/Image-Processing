% ECE5367 - P03: Image Segmentation and Detection Project
% Authors:  Joshua Williams
%           Christian Pennington
%
%
% Compatability: MATLAB release R2015a or later
%------------------------------------------------------------------------------------------------

close all; clear; clc;


%% load external matrices (for number & suit recognition later)

%numbers
load('numbers.mat','num1', 'num2', 'num3', 'num4', 'num5', 'num6', 'num7','num8','num9','num10','num11','num12','num13');
numMasks = cat(3, num1, num2, num3, num4, num5, num6, num7, num8, num9, num10, num11, num12, num13);
numNames = ["Ace", "2", "3", "4", "5", "6", "7", "8", "9", "10", "Jack", "Queen", "King"];

%suits
load('suits.mat', 'suitHeart', 'suitDiamond', 'suitClub', 'suitSpade');
suitMasks = cat(3, suitHeart, suitDiamond, suitClub, suitSpade);
suitNames = ["Hearts", "Diamonds", "Clubs", "Spades"];

%% main code
fprintf('starting...\n');

cam = openWebcam(webcamlist);


atext = text(1200, -100, 'Reading Card Number...');
btext = text(500, -100, 'Reading Card Suit...');
lastNum = '0';
lastSuit = '0';
numCount = 0;
suitCount = 0;  

while(1)
    %take a photo every 1 second
    pause(0.25);
    
    orig_img = snapshot(cam);   % take photo
    gray_img = rgb2gray(orig_img);
    
    %crop and rotate card
    justCard = straightenCrop(gray_img);    %mostly project 2 stuff
    %imshow(justCard);
    
    %subplot(2,2,1); imshow(justCard); title('Straightened Card');
    
    %do object detection on card face
    [allObjects, objectProps] = detectObjs(justCard);           %filter and get objects on card face
    if (size(objectProps, 1) > 1)
        [num_img, suit_img] = findNumSuit(allObjects, objectProps); %find objects for the number and suit
    
        %subplot(2,2,2); imshow(num_img); title('Number Seen on the Card');
        %subplot(2,2,3); imshow(suit_img); title('Suit Seen on the Card');
        
        %determine number
        index = findClosestMatch(num_img, numMasks);%returns which element in "numMasks" is the closest match to "num_img"
        numFinal = numNames(index);
        
        %determine suit
        index = findClosestMatch(suit_img, suitMasks);
        suitFinal = suitNames(index);
    
        %filter to see if our results are reliable
        if numFinal == lastNum
            numCount = numCount + 1;
        else
            numCount = 0;
        end

        if suitFinal == lastSuit
            suitCount = suitCount + 1;
        else
            suitCount = 0;
        end

        if numCount > 1
            delete(atext);
            atext = text(1200,-100, sprintf('Number: %s', numFinal));
            numCount = 0;
        end

        if suitCount > 1
            delete(btext);
            btext = text(500,-100, sprintf('Suit: %s', suitFinal));
            suitCount = 0;
        end

        %fprintf('\nGuess: Card is a %s of %s?\n', numFinal, suitFinal);
        %fprintf('numCount: %i, suitCount: %i\n', numCount, suitCount);


        lastNum = numFinal;
        lastSuit = suitFinal;
    end
end

clear; close all;

%% save (when "training" the program)
% suitSpade = suit_img;
% fprintf('Saving suit matrix...\n');
% save('suits.mat','suitSpade', '-append');


%% functions

function justCard = straightenCrop(orig_img)
    %blur
    h = ones(3,3)./9;
    img = imfilter(orig_img,h);

    %up the contrast
    img = imadjust(img);

    %quantize
    tmp = img;
    thresh = 130;
    tmp(img < thresh) = 0;
    numLevels = 5;
    delta = round( (255-thresh)./numLevels );
    while thresh < 255
        tmp(img > thresh) = thresh;
        thresh = thresh + delta;
    end
    img = tmp;

    
    %filter out noise in image
%     img = wiener2(img,[3 3]);
 
    % Median Filter the image:
%     img = medfilt2(img, [3 3]);

    
    % determine second derivative using convolution
    h = [0 1 0; 1 -4 1; 0 1 0];         %"Laplacian Kernel"
    img = edge(img,'zerocross',h);      % convolve and find zerocrossings

    img = bwareaopen(img,15); %remove small objects
    img = imfill(img, 'holes'); % fill in holes
    img = bwareaopen(img,25); %remove small objects

    
    % identify separate objects
    labeledImage = bwlabel(img, 8); %makes all pixels in a blob the same value (but unique to other blobs)
    blobProps = regionprops(labeledImage, img, 'all'); %extract information calculated during bwlabel()

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
    % Median Filter the image:
    img = medfilt2(justCard, [3 3]);

    %up the contrast
    img = imadjust(img);

    %select darker pixels only
    img = img < 150;
    img = bwareaopen(img,15); %remove small objects
    img = imfill(img, 'holes'); % fill in holes
    
    
    % identify separate objects
    labeledImage = bwlabel(img, 8); %makes all pixels in a blob the same value (but unique to other blobs)
    blobProps = regionprops(labeledImage, img, 'all'); %extract information calculated during bwlabel()
    %numberOfBlobs = size(blobProps, 1);
    
    
    % filter out objects that are super long shaped
    allBlobEccentricities = [blobProps.Eccentricity];
    allowableIndexes = (allBlobEccentricities < 0.87); % logical list... [0, 1, 1, 0]
    keeperIndexes = find(allowableIndexes); % get indices of keeper blobs

    % Extract the blobs that meet our criteria (Result will be an image)
    keeperBlobsImage = ismember(labeledImage, keeperIndexes);

    % Re-label, but only the keeper blobs
    labeledImage = bwlabel(keeperBlobsImage, 8);
    blobProps = regionprops(labeledImage, keeperBlobsImage, 'all');
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

function cam = openWebcam(webcamlist)
    %get list of webcams
    camlist = webcamlist;
    fprintf('found %i webcam(s), ', length(camlist));
    fprintf('connecting to %s\n', camlist{1});

    %open webcam preview in figure window
    cam = webcam(1);
    fig = figure();
    fig.Name = 'Webcam Preview';
    
    %make sure image isn't warped
    ax = axes(fig);         %get axes from figure
    frame = snapshot(cam);  %take a temp snapshot to what size figure window should be
    camout = image(ax,zeros(size(frame),'uint8')); 
    axis(ax, 'image');

    fprintf('opening webcam preview...\n');
    previewImageObject = preview(cam, camout);

    %Set image to be a mirror image view
    previewImageObject.Parent.XDir='reverse';

    setappdata(fig, 'cam', cam);
end