function [f] = Scint_City_fun(lambda,theta,d,n,i_scint,dz_Nz, coupled, is_Gz)
%% Constants and Variables

num_layers = length(n);

if length(lambda) == 1
    N_3rdD = length(theta);
elseif length(theta) == 1
    N_3rdD = length(lambda);
else
    fprintf('theta and lambda can`t be both vectors\n');
    f = [];
    return;
end

d_tot   = [0, d, 0];

if is_Gz
    max_Nz = dz_Nz;
    b = zeros(length(i_scint), max_Nz); % distance of dipole from closer bottom interface
    for i = 1:length(i_scint)
        delta = d(i_scint(i) - 1) / max_Nz;
        b(i, 1:max_Nz) = linspace(delta, d(i_scint(i) - 1) - delta, max_Nz);
    end
else
    dz = dz_Nz;
    Nz = zeros(1, length(i_scint));
    max_Nz = ceil(max(d(i_scint - 1)) / dz);
    b = zeros(length(i_scint), max_Nz); % distance of dipole from closer bottom interface
    for i = 1:length(i_scint)
        Nz(i) = length(dz : dz : d(i_scint(i) - 1) - dz);
        b(i, 1:Nz(i)) = dz : dz : d(i_scint(i) - 1) - dz;
    end
end

b = repmat(b, 1, 1, N_3rdD);
d_scint = repmat(d_tot(i_scint).', 1, max_Nz, N_3rdD);
n_scint = repmat(n(i_scint).',     1, max_Nz, N_3rdD);
top = d_scint - b; % distance of dipole from closer top interface
mask = b ~= 0;

%% Effective Fresnel's coeff calculation
% R_eff and T_eff return a vector

last_k = 2*pi * n(end) ./ lambda;
u      = last_k * sin(theta);

% R_eff(eps,d,u,lambda,type)  
r.s_up   = R_eff(n(2:end).^2,[d(2:end) 0],u,lambda,'s');
r.p_up   = R_eff(n(2:end).^2,[d(2:end) 0],u,lambda,'p');
r.s_down = R_eff(fliplr(n(1:end - 1).^2),fliplr([0 d(1:end - 1)]),u,lambda,'s');
r.p_down = R_eff(fliplr(n(1:end - 1).^2),fliplr([0 d(1:end - 1)]),u,lambda,'p');

% T_eff(eps,d,u,r,lambda,type)  
t.s_up   = T_eff(n(2:end).^2,[d(2:end) 0],u,r.s_up,lambda,'s');
t.p_up   = T_eff(n(2:end).^2,[d(2:end) 0],u,r.p_up,lambda,'p');

% Choosing only the relevant coeffs
r.s_up   = r.s_up(i_scint - 1, :);            r.p_up   = r.p_up(i_scint - 1, :); 
r.s_down = r.s_down(num_layers - i_scint, :); r.p_down = r.p_down(num_layers - i_scint, :);
t.s_up   = t.s_up(i_scint - 1, :);            t.p_up   = t.p_up(i_scint - 1, :);

r.s_up   = permute(repmat(r.s_up,   1, 1, max_Nz), [1 3 2]); 
r.p_up   = permute(repmat(r.p_up,   1, 1, max_Nz), [1 3 2]); 
r.s_down = permute(repmat(r.s_down, 1, 1, max_Nz), [1 3 2]); 
r.p_down = permute(repmat(r.p_down, 1, 1, max_Nz), [1 3 2]); 
t.s_up   = permute(repmat(t.s_up,   1, 1, max_Nz), [1 3 2]); 
t.p_up   = permute(repmat(t.p_up,   1, 1, max_Nz), [1 3 2]); 

%%  Calculating the enhancment for all scintillator's layers

if length(lambda) == 1
    last_k  = repmat(2*pi * n(end)       / lambda, length(i_scint), max_Nz, N_3rdD);
    k_scint = repmat(2*pi * n(i_scint).' / lambda, 1, max_Nz , N_3rdD);
else
    lambda_mat = permute(repmat(lambda, length(i_scint), 1 , max_Nz), [1 3 2]);
    last_k  = repmat(2*pi * n(end), length(i_scint), max_Nz, N_3rdD)./lambda_mat;
    k_scint = repmat(2*pi * n(i_scint).', 1, max_Nz , N_3rdD)./lambda_mat;
end


u = repelem(u, length(i_scint) * max_Nz);
u = reshape(u, length(i_scint), max_Nz, N_3rdD);

last_l  = sqrt(last_k.^2  - u.^2);
l_scint = sqrt(k_scint.^2 - u.^2);

if coupled

    T_par_s_up =  t.s_up .* (1 + r.s_down .* exp(2i * l_scint .* b)) ./ (1 - r.s_down .* r.s_up .* exp(2i * l_scint .* d_scint)) ;
    T_par_p_up =  t.p_up .* (1 + r.p_down .* exp(2i * l_scint .* b)) ./ (1 - r.p_down .* r.p_up .* exp(2i * l_scint .* d_scint)) ;
    T_perp_p_up = t.p_up .* (1 - r.p_down .* exp(2i * l_scint .* b)) ./ (1 - r.p_down .* r.p_up .* exp(2i * l_scint .* d_scint)) ;
    temp = (1/4) * (last_k .* last_l.^2 ./ k_scint.^3);

    Par_s_up  = temp .* ( k_scint.^2 .* (n_scint / n(end)).^2 .* abs( T_par_s_up ./ l_scint ).^2) ;
    Par_p_up  = temp .* abs( T_par_p_up ).^2 ;
    Perp_p_up = temp .* u.^2 .*abs( T_perp_p_up ./ l_scint ).^2 ;

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

f = squeeze((sum(sum(Gz(f,d_tot(i_scint).',is_Gz), 2), 1)));

end

function G = Gz(f,d,is_Gz)
    if is_Gz
        G = permute(f, [3 2 1]);
        for i = 1:length(d)
            G(:,:,i) = G(:,:,i) * d(i);
        end
        G = permute(G, [3 2 1]);
    else
        G = f;
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
        r = zeros(max(length(lambda), length(u)));
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
        r = [(r12 + R(1,:) .* exp(2i*l2*d(1))) ./ (1 + r12 .* R(1,:) .* exp(2i*l2*d(1))); R];
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
        t = ones(max(length(lambda), length(u)));
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
            t12 = 2*l1*sqrt(eps(1)*eps(2))./(l2*eps(1) + l1*eps(2));% TM
        else
            t = [];
            return
        end

        T = T_eff(eps(2:end),d(2:end),u,R(2:end,:),lambda,type);
        t = [(t12 .* T(1,:) .* exp(1i*l2*d(1))) ./ (1 + r12 .* R(2,:) .* exp(2i*l2*d(1))); T];
    end
end
