%ECE 5367: Image Processing
%Project 1
%Author: Joshua Williams

close all; clear; clc;

%path = 'C:\Users\joshu\Documents\GitHub\ImageProcessing\P01\ImageSet1\';
path = input('Enter file path: ', 's');

images = dir(fullfile(path,'*.jpg')); % array of all .jpg images in path
numImages = length(images);

disp('Displaying images, press any key to move to next image');

for k = 1:numImages
    
    % read RGB image
    RGBimg = imread(fullfile(path,images(k).name));
    
    % convert to HSV
    HSVimg = rgb2hsv(RGBimg);
    HSVsat = HSVimg(:, :, 2); %saturation
    
    % find mean saturation value
    meanS = mean(HSVsat(:));
    
    % determine if it is day or night based on the mean saturation value
    if (meanS < 0.1)    %0.1 threshold determined experimentally
        category = 'night';
    else
        category = 'day';
    end
    
    imshow(RGBimg);
    title(category);
    
    pause();
end

disp('done!');