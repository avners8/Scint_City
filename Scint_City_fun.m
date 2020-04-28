function [f] = Scint_City_fun(lambda,theta,d,n,i_scint,coupled,control)
% This function calculate the emission rate enhancment of a given
% structure.
% Inputs:
%     theta   - the solid angle in which we compute the emission rate (a scalar or a vector)
%     lambda  - the wavelength in which we compute the emission rate (a scalar or a vector)
%     d       - a vector of the thicknesses of the structure
%     n       - a vector of the refractive indices of the structure
%     i_scint - a vector of the indices of the scintillator layers
%     coupled - if true, compute the emission rate of the coupled wave(outside the structure) and the wave inside the device.
%     control - a struct that contains computation parameters:
%                   1) is_Gz - There are two mode of computation - Gz or dz. is_Gz choose between those modes.
%                      Gz - distributes a fix number of emitters in each layer of scintillator, then weight the contribution of each dipole to the total emission, with the Gz function.
%                      dz - in this mode we assume a distance of dz between two emitters. Thus, thicker layers will have more emitters.
%                   2) dz_Nz - In Gz mode, this parameter is the fix number of emitters in each layer. In dz mode, this number is the distance between two emitters.
%                   3) sum_on_z - if True, the function sum all the
%                   contribution of the emitters and return the total emission rate. Otherwise, the function return the emission rate of each dipole separetly.
% Outputs:
%   f - the emission rate inside or outside the device.
%   if sum_on_z is true the dimensions of f are [N_theta X N_lambda],
%   otherwise the dimensions of f are [Scintilator layers X emitters in each layer X N_theta X N_lambda]
%% Constants and Variables

dz_Nz = control.dz_Nz; is_Gz = control.is_Gz; sum_on_z = control.sum_on_z;

num_layers = length(n);

% theta and lambda can be a single value or a vector of value
% the output is a matrix of the emission rate per angle and wavelength
N_theta  = length(theta);
N_lambda = length(lambda);

d_tot   = [0, d, 0];

% b is the distance of each emitter to the nearest bottom of the interface.
% compute_b compute this matrix 
[b, max_Nz]	= compute_b(d,i_scint,N_theta,N_lambda,is_Gz,dz_Nz);
d_scint = repmat(d_tot(i_scint).', 1, max_Nz, N_theta, N_lambda); % the vector of the thicknesses of the scintillator layers
n_scint = repmat(n(i_scint).',     1, max_Nz, N_theta, N_lambda); % the vector of the refractive indices of the scintillator layers
top = d_scint - b;                                                % distance of dipole from nearest top interface

% In dz mode, each layer can have a different layer, thus a different
% number of emitters in each layer. We assigned a matrix assuming all layers have
% the thickness of the thickest layer, then we zero the irrelevant emitters using a mask.
mask = b ~= 0;

%% Effective Fresnel's coeff calculation

% R_eff and T_eff return a vector
theta_mat  = repmat(theta, 1, 1, N_lambda);                  
lambda_mat = permute(repmat(lambda, N_theta, 1, 1), [3,1,2]);
last_k     = repmat(2*pi * n(end), 1, N_theta, N_lambda)./lambda_mat;
u          = last_k .* sin(theta_mat);

% R_eff(eps,d,u,lambda,type)  
r.s_up   = R_eff(n(2:end).^2,[d(2:end) 0],u,lambda_mat,'s');
r.p_up   = R_eff(n(2:end).^2,[d(2:end) 0],u,lambda_mat,'p');
r.s_down = R_eff(fliplr(n(1:end - 1).^2),fliplr([0 d(1:end - 1)]),u,lambda_mat,'s');
r.p_down = R_eff(fliplr(n(1:end - 1).^2),fliplr([0 d(1:end - 1)]),u,lambda_mat,'p');

% T_eff(eps,d,u,r,lambda,type)  
t.s_up   = T_eff(n(2:end).^2,[d(2:end) 0],u,r.s_up,lambda_mat,'s');
t.p_up   = T_eff(n(2:end).^2,[d(2:end) 0],u,r.p_up,lambda_mat,'p');

% Choosing only the relevant coeffs
r.s_up   = r.s_up(i_scint - 1, :, :);            r.p_up   = r.p_up(i_scint - 1, :, :); 
r.s_down = r.s_down(num_layers - i_scint, :, :); r.p_down = r.p_down(num_layers - i_scint, :, :);
t.s_up   = t.s_up(i_scint - 1, :, :);            t.p_up   = t.p_up(i_scint - 1, :, :);

r.s_up   = permute(repmat(r.s_up,   1, 1, 1, max_Nz), [1 4 2 3]); 
r.p_up   = permute(repmat(r.p_up,   1, 1, 1, max_Nz), [1 4 2 3]); 
r.s_down = permute(repmat(r.s_down, 1, 1, 1, max_Nz), [1 4 2 3]); 
r.p_down = permute(repmat(r.p_down, 1, 1, 1, max_Nz), [1 4 2 3]); 
t.s_up   = permute(repmat(t.s_up,   1, 1, 1, max_Nz), [1 4 2 3]); 
t.p_up   = permute(repmat(t.p_up,   1, 1, 1, max_Nz), [1 4 2 3]); 

%%  Calculating the enhancment for all scintillator's layers

theta_mat  = permute(repmat(theta' , 1, length(i_scint), max_Nz, N_lambda), [2 3 1 4]);
lambda_mat = permute(repmat(lambda', 1, length(i_scint), max_Nz, N_theta) , [2 3 4 1]);
k          = repmat(2*pi * n(end), length(i_scint), max_Nz, N_theta, N_lambda)./lambda_mat;
k_scint    = repmat(2*pi * n(i_scint).', 1, max_Nz , N_theta, N_lambda)./lambda_mat;
if coupled
    u      = k .* sin(theta_mat);
else % inside
    u      = k_scint .* sin(theta_mat);
end

last_l  = sqrt(k.^2  - u.^2);
l_scint = sqrt(k_scint.^2 - u.^2);

if coupled

    T_par_s_up =  t.s_up .* (1 + r.s_down .* exp(2i * l_scint .* b)) ./ (1 - r.s_down .* r.s_up .* exp(2i * l_scint .* d_scint)) ;
    T_par_p_up =  t.p_up .* (1 + r.p_down .* exp(2i * l_scint .* b)) ./ (1 - r.p_down .* r.p_up .* exp(2i * l_scint .* d_scint)) ;
    T_perp_p_up = t.p_up .* (1 - r.p_down .* exp(2i * l_scint .* b)) ./ (1 - r.p_down .* r.p_up .* exp(2i * l_scint .* d_scint)) ;
    temp = (1/4) * (k .* last_l.^2 ./ k_scint.^3);

    Par_s_up  = temp .* ( k_scint.^2 .* (n_scint / n(end)).^2 .* abs( T_par_s_up ./ l_scint .* exp(1i * l_scint .* top)).^2) ;
    Par_p_up  = temp .* abs( T_par_p_up .* exp(1i * l_scint .* top)).^2 ;
    Perp_p_up = temp .* u.^2 .*abs( T_perp_p_up ./ l_scint .* exp(1i * l_scint .* top)).^2 ;

    f = (Par_s_up + Par_p_up + Perp_p_up) .* mask;

else % Purcel Factor

    R_perp_p_up   = (r.p_down .* (r.p_up   .* exp(2*1i.*l_scint.*top) - 1)) ./ (1 - r.p_down .* r.p_up .* exp(2*1i.*l_scint.* d_scint ));
    R_perp_p_down = (r.p_up   .* (r.p_down .* exp(2*1i.*l_scint.*b)   - 1)) ./ (1 - r.p_down .* r.p_up .* exp(2*1i.*l_scint.* d_scint ));
    R_par_s_up    = (r.s_down .* (r.s_up   .* exp(2*1i.*l_scint.*top) + 1)) ./ (1 - r.s_down .* r.s_up .* exp(2*1i.*l_scint.* d_scint ));
    R_par_s_down  = (r.s_up   .* (r.s_down .* exp(2*1i.*l_scint.*b)   + 1)) ./ (1 - r.s_down .* r.s_up .* exp(2*1i.*l_scint.* d_scint ));
    R_par_p_up    = (r.p_down .* (r.p_up   .* exp(2*1i.*l_scint.*top) + 1)) ./ (1 - r.p_down .* r.p_up .* exp(2*1i.*l_scint.* d_scint ));
    R_par_p_down  = (r.p_up   .* (r.p_down .* exp(2*1i.*l_scint.*b)   + 1)) ./ (1 - r.p_down .* r.p_up .* exp(2*1i.*l_scint.* d_scint ));
    temp = 1 ./ (2 * k_scint .* l_scint);

    Par_s_up  = temp .* u.^2        .* (1 + R_perp_p_up.*exp(2*1i.*l_scint.*b) + R_perp_p_down.*exp(2*1i.*l_scint.*top) );
    Par_p_up  = temp .* k_scint.^2  .* (1 + R_par_s_up .*exp(2*1i.*l_scint.*b) + R_par_s_down .*exp(2*1i.*l_scint.*top) );
    Perp_p_up = temp .* l_scint.^2  .* (1 + R_par_p_up .*exp(2*1i.*l_scint.*b) + R_par_p_down .*exp(2*1i.*l_scint.*top) );

    f = real(Par_s_up + Par_p_up + Perp_p_up) .* mask;

end

if sum_on_z
    f = squeeze((sum(sum(Gz(f,d_tot(i_scint).',is_Gz), 2), 1)));
end

end

function [r] = R_eff(eps,d,u,lambda,type)
% R_eff Calculates the special return coefficient
% 
% Input
%   eps - Vector of the premitivity coefficients.
%   d - Vector of the layers widths.
%   u - Planer wave vector size.
%   lambda - Wavelength.
%   type - s (TE) or p (TM).
%
% Returns
%   r = The return coefficients.
%
% Drawing
%             d1      d2      d3       d4     d5
%           ------- ------- ------- ------- -------     ...
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%    R*--> |       |       |       |       |       |    ...
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%    eps1    eps2     eps3    eps4    eps5    eps6      ...


    if (length(eps)==1)
        r = zeros(size(u));
    else
        k1 = 2*pi*(1./lambda)*sqrt(eps(1));
        k2 = 2*pi*(1./lambda)*sqrt(eps(2));

        l1 = sqrt(k1.^2-u.^2); % kz1
        l2 = sqrt(k2.^2-u.^2); % kz2

        if (type == 's')
            r12 = (l1 - l2)./(l1 + l2); % TE
        elseif (type == 'p')
            r12 = (l2*eps(1) - l1*eps(2))./(l2*eps(1) + l1*eps(2)); % TM
        else
            r = [];
            return
        end
        R = R_eff(eps(2:end),d(2:end),u,lambda,type);
        r = [(r12 + R(1,:,:,:) .* exp(2i*l2*d(1))) ./ (1 + r12 .* R(1,:,:,:) .* exp(2i*l2*d(1))); R];
    end
end

function [t] = T_eff(eps,d,u,R,lambda,type)
% R_eff Calculates the special return coefficient
% 
% Input
%   eps - Vector of the premitivity coefficients.
%   d - Vector of the layers widths.
%   u - Planer wave vector size.
%   lambda - Wavelength.
%   type - s (TE) or p (TM).
%
% Returns
%   r = The return coefficient.
%
% Drawing
%             d1      d2      d3       d4     d5
%           ------- ------- ------- ------- -------     ...
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%    R*--> |       |       |       |       |       |    ...
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%    eps1    eps2     eps3    eps4    eps5    eps6      ...


    if (length(eps)==1)
        t = ones(size(u));
    else
        k1 = 2*pi*(1./lambda)*sqrt(eps(1));
        k2 = 2*pi*(1./lambda)*sqrt(eps(2));

        l1 = sqrt(k1.^2-u.^2); % kz1
        l2 = sqrt(k2.^2-u.^2); % kz2

        if (type == 's')
            r12 = (l1 - l2)./(l1 + l2); % TE
            t12 = 2*l1./(l1 + l2);      % TE
        elseif (type == 'p')
            r12 = (l2*eps(1) - l1*eps(2)) ./(l2*eps(1) + l1*eps(2));% TM
            t12 = 2*l1.*sqrt(eps(1)*eps(2))./(l2*eps(1) + l1*eps(2));% TM
        else
            t = [];
            return
        end

        T = T_eff(eps(2:end),d(2:end),u,R(2:end,:,:,:),lambda,type);
        t = [(t12 .* T(1,:,:,:) .* exp(1i*l2*d(1))) ./ (1 + r12 .* R(2,:,:,:) .* exp(2i*l2*d(1))); T];
    end
end
