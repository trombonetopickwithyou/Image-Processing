%Image Processing Project 4
%Joshua Williams

close all; clear; clc;


%% load image
filename = 'Proj4.tif'; %program looks in project directory
img = im2double(imread(filename));

subplot(2,3,1);
imshow(img);
title('original image');
hold on;

%% fourier transform image
%F = log(1+abs(fftshift(fft2(img))));
F = fft2(img);
F = fftshift(F);
F = log(1 + F);


subplot(2,3,4);
surf(abs(F));
title('fourier transform');

%% construct filter (ellipse in freq domain)
[numrows, numcols] = size(F);
centerRow = round(numrows/2);
centerCol = round(numcols/2);

%cutoff "frequencies" in the x and y directions
cutoffY = 25;
cutoffX = 10;


% create an ellipse filter in freq domain
theta = linspace(0,2*pi);
col = centerCol + cutoffX.*cos(theta);
row = centerRow + cutoffY.*sin(theta);

H = zeros(numrows, numcols);
for i=1:numrows
    for j = 1:numcols
        if(inpolygon(i,j,row,col))  %make a filled in white oval
            H(i,j) = 1;
        end
    end
end

subplot(2,3,2);
imshow(H);
title('filtering mask to be used');



% apply filter
F = F.*H;
subplot(2,3,5);
surf(abs(F));
title('filtered fourier transform');

%% inverse fft
%out = abs(ifft2(F));
F = exp(F) - 1;
F = ifftshift(F);
out = real(ifft2(F));


%% filter image noise (in spatial domain)
%out = imadjust(out, [0.45 0.8]);
out = imadjust(out, [0.25 0.65]);    %adjust contrast ([0-1 lwr bound, 0-1 upr bound])

out = medfilt2(out, [5 5]);                         %median filter
out = imsharpen(out, 'Radius', 5, 'Amount', 2);   %unsharp masking


subplot(2,3,3);
imshow(out);
title('output image');


