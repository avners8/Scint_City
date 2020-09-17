function G = Gz(f,d_tot,i_scint,control)
    if control.is_Gz
        if control.is_exponential_absorption
            G = permute(f, [3 4 1 2]);
            d = d_tot(2:(end-1));
            bottom_distance = compute_b(d,i_scint,1,1,control.is_Gz,control.dz_Nz);
            if control.plot_exponential_absorption 
                figure();
            end
            for i = 1:length(i_scint) % Iterating over all scintillator layers
                for j = 1:control.dz_Nz % Iterating over all emitters in a layer
                    full_layers_scint_distance = sum( d_tot(i_scint( 1:(i-1) ) )); 
                    scint_distance = full_layers_scint_distance + bottom_distance(i,j,1,1);
                    other_distance = sum( d_tot( 1: (i_scint(i) - 1) )) - full_layers_scint_distance;
                    G(:,:,i,j) = G(:,:,i,j) * exp(- scint_distance/control.absorption_scint - other_distance/control.absorption_other);
                    if control.plot_exponential_absorption 
                        plot(scint_distance + other_distance, exp(- scint_distance/control.absorption_scint - other_distance/control.absorption_other), 'o'); hold on;
                    end
                end
            end
            if control.plot_exponential_absorption 
                graphParams('Exponential absorption profile in the structure', '$Thickness\ [nm]$', '$Weight$');
            end
            G = permute(G, [3 4 1 2]);
        else % uniform distribution of emitters
            d_scint = d_tot(i_scint);
            G = permute(f, [4 2 3 1]);
            for i = 1:length(d_scint)
                G(:,:,:,i) = G(:,:,:,i) * d_scint(i);
            end
            G = permute(G, [4 2 3 1]);
        end
    else
        G = f;
    end
end

function graphParams(ptitle, pxlabel, pylabel) 
    grid on;
    title(ptitle); xlabel(pxlabel); ylabel(pylabel);
    set(gca, 'FontSize', 14); set(gcf,'color','w'); set(gca,'linewidth',2.5);
end