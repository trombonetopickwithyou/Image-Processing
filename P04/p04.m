%Image Processing Project 4
%Joshua Williams

close all; clear; clc;


%% load image
filename = 'Proj4.tif'; %program looks in project directory
img = im2double(imread(filename));

subplot(2,2,1);
imshow(img);
title('original image');
hold on;

%% fourier transform image
%F = log(1+abs(fftshift(fft2(img))));
F = fft2(img);
F = fftshift(F);
F = log(1 + F);

subplot(2,2,2);
surf(abs(F));
title('fourier transform');

%% construct filter
cutoffmin = 0;
cutoffmax = 30;

[numrows, numcols] = size(F);
centerRow = round(numrows/2);
centerCol = round(numcols/2);

H = zeros(numrows, numcols);

for i=1:numrows
    for j = 1:numcols
        distance = dist(i,j,centerRow,centerCol);
        
        if(distance < cutoffmin)    
            H(i,j) = 0.1;
        elseif(distance < cutoffmax)
            H(i,j) = 1;            
        else
            H(i,j) = 0.5;%max(1-0.02*(distance-cutoffmax),0);
        end
    end
end
H(centerRow, centerCol) = 0;

F = F.*H;
subplot(2,2,4);
surf(abs(F));
title('filtered fourier transform');

%% extract periodic signal from image (filter image)
%out = abs(ifft2(F));
F = exp(F) - 1;
F = ifftshift(F);
out = abs(ifft2(F));

out = imadjust(out, [0.2 0.6]);
out = imsharpen(out, 'Radius', 0.2, 'Amount', 2);

% minval = min(out(:));
% maxval = max(out(:));
% out = 255.*(out - minval)/(maxval - minval);

subplot(2,2,3);
imshow(out);
title('output image');



%% functions
function out = dist(x1,y1,x2,y2)
    out = sqrt((x1-x2).^2 + (y1-y2).^2);
end


