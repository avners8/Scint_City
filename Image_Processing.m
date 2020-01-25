function [p_im,mse,mtf] = Image_Processing(im,psf,total_size,ImageProcessing_Params)
% This function model the 1D PhC as a filter, and process the image with it
% Input arguments:
%   im    - image to be processed
%   psf   - the psf (calculated by Scint-City)
%   k     - the reduction parameter of the psf
%   total_size - total size of the structure
% Output arguments:
%   p_im  - processed image
%   mse   - mse between original image and processed image
%   mtf   - the MTF (FT of psf)
%
%   p_im = (im * N_e) * H + N_g

    p_im = im;
    
    DR = ImageProcessing_Params.DR;
    nbins = ImageProcessing_Params.nbins;
    sigma_onion_noise = ImageProcessing_Params.sigma_onion_noise;
    
    % Calculating the MTF (in logarithmic scale)
    mtf = log(1 + abs(fftshift(fft(psf))));
    
    % Adding onion noise
    onion_filter = normpdf(((-nbins/2):(nbins/2)), 0, sigma_onion_noise);
    psf = conv(psf, onion_filter,'same'); 
    
    psf_2D  = conv2(psf.', psf); 
    psf_2D  = DR * psf_2D  / sum(sum(psf_2D) );
    
    % Adding inefficiency noise
    mu = 3 * 1e6; 
    p  = tanh(mu * total_size); % the probability that a pixel is not detected
    bernuli_mat = (rand(size(p_im)) < p)*1;
    p_im = double(p_im) .* bernuli_mat;
    
    % Blurring the image
    p_im = imfilter(p_im, psf_2D,  'replicate');
    
    % Adding geussian noise
    sigma_noise = 0.5;
    noise_g     = sigma_noise * randn(size(p_im));
    p_im = p_im + noise_g;
    
    p_im = uint8(p_im);
    
    % Calculating the mse
    mse = immse(im, p_im);

end
