function [ F ] = fft_vis(I)

    % Input:
    % I: the input image
    % Output:
    % F: 2D Fourier transform in a form amenable to visualization
    
    % Please follow the instructions in the comments to fill in the missing commands

    % 1) Apply Fourier transform to the image (MATLAB command fft2)
A = fft2(I);
    % 2) Shift the spectrum (MATLAB command fftshift)
B = fftshift(A);
    % 3) Take the absolute value
C = abs(B);
    % 4) Add 1 and take the log (for visualization)
F = log(C+ 1);
end