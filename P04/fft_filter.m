% Matlab Program to demonstrate the "Low pass Filtering of an image using
% 2D-DFT"
clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm for filtering in the frequency Domain
% Step1: Given ain input image f(x,y) of size M x N, obtain the pading
% parameters P and Q. Typically, we select P = 2M and Q = 2N
% Step2: Form a padded image fp(x,y) of size P X Q by appending the
% necessary number of zeros to f(x,y).
% Step3: Multiply fp(x,y) by (-1)^(x+y)
% Step4: Compute the DFT, F(u,v) of the image from Step 3
% Step5: Generate a Real, Symmetric Filter Function H(u,v) of size P X Q
% with center at coordinates (P/2,Q/2), 
% Step 6:Form the product G(u,v) = H(u,v)F(u,v) using array multiplication
% Obtain the processed image 
% Step 7: gp(x,y) = {real{inverse DFT[G(u,v)]}(-1)^(x+y)
% Step 8: Obtain the final processed result g(x,y) by extracting the M X N region
% from the top, left quadrant of gp(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorith Implementation

%% read in image
img = im2double(imread('Proj4.tif'));
[m,n] = size(img);

imshow(img);title('original image');

%% compute fourier transform
% Computing the 2D DFT using "fft2" matlab command
%F = abs(fft2(img));
F = log(1+abs(fftshift(fft2(img))));
figure;imshow(F);title('2D DFT of the pre processed image');

%% generate filter
[x,y] = freqspace([m n],'meshgrid');
z = zeros(m,n);

%calculate distance from center at each point
for i = 1:m
    for j = 1:n
        z(i,j) = sqrt(x(i,j).^2 + y(i,j).^2);
    end
end

% Choosing the Cut off Frequency (generating circle mask)
H = zeros(m,n);
for i = 1:m
    for j = 1:n
        if z(i,j) <= 0.4  % here 0.4 is the cut-off frequency of the LPF
            H(i,j) = 1;
        else
            H(i,j) = 0;
        end
    end
end
figure;imshow(H);title('Low Pass Filter Mask');

%% filter the image
% e : the 2D DFT output of pre processed image
% H : the mask for Low Pass Filter

h1 = F.*H;
figure;imshow(h1);title('Low passed output');

%% get inverse fft
%gp(x,y) = {real{inverse DFT[G(u,v)]}(-1)^(x+y)
out = abs(ifft2(h1));

figure;imshow([img out]);title('input image                 output image');