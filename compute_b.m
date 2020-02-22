function [b, max_Nz] = compute_b(d,i_scint,N_theta,N_lambda,is_Gz,dz_Nz)
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

    b = repmat(b, 1, 1, N_theta, N_lambda);
end