function [x,psf] = psf_computation(d,n,i_scint,lambda,control,ImageProcessing_Params)
% This function compute the psf of a structure
% Input arguments:
%   d           - the thicknesses of the structure
%   n           - the refractive indices
%   i_scint     - indices of the scintillating layer
%   control     - general control
%   IP_Params   - General parameter we need to set in our model
% Output arguments:
%   x           - the array of pixel (x = Ts(theta))
%   psf         - psf of the structure (a 1D vector)

    N_theta = 10001;
    theta = linspace(-pi/2, pi/2, N_theta);
    control.sum_on_z = false;
    coupled = true;
    is_Gz = control.is_Gz;
    d_tot   = [0, d, 0];
    
    h = ImageProcessing_Params.h;
    max_distance = ImageProcessing_Params.max_distance;
    
    % Computing f
    f = Scint_City_fun(lambda,theta,d,n,i_scint,coupled,control);
    
    % calculating x
    [z, max_Nz] = compute_b(d,i_scint,N_theta,1,control.is_Gz,control.dz_Nz);
    x      = zeros(size(f));
    mask   = zeros(size(x));
    mask_h = zeros(size(x));
    for i = 1 : length(i_scint)
       for b_i = 1 : max_Nz
           b = z(i, b_i,1);
           theta1 = asin((n(i) ./ n(end)) .* sin(theta));
           mask(i, b_i, :) = (theta1 == real(theta1));           
           [x(i, b_i, :),mask_h(i, b_i, :)] = Ts(d,n,i_scint(i),i_scint(i),b,h,theta1 .* (squeeze(mask(i, b_i, :))).');
       end
   end
    
   right_limit = max_distance;
   left_limit = -right_limit;
   nbins = ImageProcessing_Params.nbins;
   quantized_x = linspace(left_limit,right_limit,nbins);
   g_f = Gz(f .* mask .* mask_h,d_tot.',i_scint,control);
   psf = zeros(size(quantized_x));
   for i = 1 : length(i_scint)
     for b_i = 1 : max_Nz
        index = quantiz(squeeze(x(i, b_i, :)), quantized_x);
           index(index == 0) = 1;
           for bin = 1:nbins
                psf(bin) = psf(bin) + sum(g_f(i, b_i, (index == bin)));
           end
     end
   end
    
    psf(1) = 0; psf(end) = 0;
    x = quantized_x;
end

function [x,mask_h] = Ts(d,n,i_scint,layer_index,z,h,theta)
    if layer_index == i_scint
        [x,mask_h] = Ts(d, n, i_scint, layer_index + 1, z, h, theta);
        mask_h = imag(theta) == 0;
        x = x + z .* tan(theta) .* mask_h;
    else
        theta2 = asin((n(layer_index-1) ./ n(layer_index)) .* sin(theta));
        if layer_index - 1 > length(d)
            mask_h = imag(theta2) == 0;
            x = h .* tan(theta2) .* mask_h;
        else
            [x,mask_h] = Ts(d, n, i_scint, layer_index + 1, z, h, theta2);
            mask_h = imag(theta2) == 0;
            x = x + d(layer_index - 1) .* tan(theta2) .* mask_h;
        end
    end
end
