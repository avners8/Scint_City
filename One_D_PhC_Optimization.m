function One_D_PhC_Optimization(                ...
                                theta,          ...
                                lambda,         ...
                                mu_lambda,      ...
                                sigma_lambda,   ...
                                pairs,          ...
                                n_substrate,    ...
                                n_scint,        ...
                                n_other,           ...
                                d0,             ...
                                n,              ...
                                i_scint,        ...
                                control,          ...
                                random_iterations, ...
                                total_size,     ...
                                lower_bound,    ...
                                upper_bound,    ...
                                constraint, ...
                                optimize,       ...
                                optimization_mode, ...
                                ImageProcessing_Params, ...
                                save_fig,       ...
                                dir_name,       ...
                                simulation_params)
                            
    % This function is the main function that is in charge of the
    % structure's optimization process.
    % In addition to the optimization, the function show the results for
    % the given strcutures.
    % Inputs:
    % 
    % theta         - the solid angle in which we optimize (a scalar or a vector)
    % lambda        - the wavelength in which we optimize (a scalar or a vector)
    % mu_lambda     - the expected value E of the spectral emission
    % sigma_lambda  - the std of the spectral emission
    % pairs         - the structure is based on a substrate, then "pairs" pairs of scintillator and other dialectric layers (a scalar or a vector).
    % n_substrate   - the refractive ines of the substrate
    % n_scint       - the refractive ines of the scintillator layers
    % n_other       - the refractive ines of the other dialetric materials layers
    % d0            - if the structure is a known structure, d0 is the vector of thicknesses.
    % n             - if the structure is a known structure, n is the vector of refractive indices.
    % i_scint       - if the structure is a known structure, i_scint is the vector of indices of the scintillator layers.
    % is_Gz         - there are two mode of computation - Gz or dz. is_Gz choose between those modes. 
    % random_iterations - the function performs the optimization process several times with random starting points to improve the results. This choose the number of itarations.
    % total_size    - Used to define the constraints on the optimization process, and to choose the bulk structure to compare to.
    % constraint    - Choose the constraint mode. Explained in the GUI
    % optimize      - if true, it perform the optimization. Otherwise, it expect a known structure.
    % ImageProcessing_Params - the parameters for the image processing - nbins, max distance of propagation ...
    % save_fig      - if true saved all the figs 
    % dir_name      - if save_figs is true, save all figs to specified directory.     
    % simulation_params

% Drawing
%             d1      d2      d3       d4     d5
%           ------- ------- ------- ------- -------     ...
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%          |       |       |       |       |       |
% substrate| scint | other | scint | other | scint |    ...
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%          |       |       |       |       |       |
%  n_subst. n_scint n_other n_scint n_other n_scint     ...

    %% Spectral distribution
    
    % if the std of the spectral distribution is 0, we calculate for only one wavelength.
    if sigma_lambda == 0
        Y_orig = zeros(size(lambda));
        Y_orig(lambda == mu_lambda) = 1;
    else
        Y_orig = normpdf(lambda, mu_lambda, sigma_lambda);
    end
    if length(theta) == 1
        Y      = Y_orig;
    else
        Y      = repmat(Y_orig, length(theta), 1);
    end
    dz_Nz = control.dz_Nz; is_Gz = control.is_Gz;  control.sum_on_z = true;

    %% General parameters

    % Saving the original values of theta and lambda. Will bu used when plotting.
    theta_orig = theta;
    lambda_orig = lambda;
    num_theta_plot = 1001;

    % Defining the colors in the plots.
    lighter = - 0.071;

    Green = [0.4660 0.6740 0.1880] - lighter;
    Red   = [0.6350 0.0780 0.1840] - lighter;
    Gold  = [0.9290 0.6940 0.1250] - lighter;

    Green_Sim = [0.4660 0.6740 0.1880] + lighter;
    Red_Sim   = [0.6350 0.0780 0.1840] + lighter;
    Gold_Sim  = [0.9290 0.6940 0.1250] + lighter;
    
    % Adding the current directory to the provided folder name.
    % Make sure the folder was created before using this function.
    dir_name = [pwd, '\', dir_name];


    %% Calculating the minimum using Matlab Optimization Toolbox

    % Scint_City compute the emission rate enhancment of the coupled wave
    % (outside the structure) and the wave inside the device.
    coupled = true;
    inside  = false;
	
    % If Y is a single wavelength, transforming Y to a scalar
    if sigma_lambda == 0
       lambda = mu_lambda;
       Y      = 1;
    end
    
    if ((optimize == 0) & control.load_data_from_dir)
        load(dir_name + "\optimized_d.mat");
        load(dir_name + "\n.mat");
        load(dir_name + "\i_scint.mat"); 
    end
    
    if optimize
        n_bulk  = [n_substrate n_scint 1];
    else
        n_bulk  = [n(1) n(i_scint(1)) n(end)];
    end
    
    % There is two modes of operation.
    % First mode is optimize: in this mode the function get the constraints
    % and the general parameters of the structure. The function find an
    % optimized structure in term of emission rate.
    % Secondly, the function get a known structure, and compute it's
    % emission rate enhancment and efficiency. 
    
    % This section is relavant only for the optimization process.
    if optimize
        y_max_max = 0;
        F_max_max = 0;
        pairs_max = pairs(1);
        
        % Iterating on number of layers in the structure
        for p = 1:length(pairs)
            
            % Creating n and i_scint, for the specific structure with
            % "pairs" pairs.
            n       = [n_substrate, repmat([n_scint, n_other], 1 , pairs(p)), 1];
            i_scint = 2:2:2 * pairs(p);
            
            i_other = zeros(1,length(n) - length(i_scint) - 2);
            count_scint = 1; count_other = 1;
            for j = 2:length(n) - 1
                if j - 1 == i_scint(count_scint)
                    count_scint = count_scint + 1;
                else
                    i_other(count_other) = j + 1;
                    count_other = count_other + 1;
                end
            end
            
            % Calculating constraints
            % starting point
            even = 0.5 * (1+(-1).^(1:(2*pairs(p))));
            odd  = 0.5 * (1+(-1).^(0:(2*pairs(p) - 1)));

            % Defining the set of linear constraints - Ad=b
            % There are 3 possible constraints.
            
            A = [];
            b = [];
            if constraint == 1      % 1: Total scint and Total other equal Total size
                Aeq = [even;odd];
                beq = [total_size,total_size];
            elseif constraint == 2  % 2: Each layer is between 0 and Total size
                Aeq = [];
                beq = [];
            elseif constraint == 3  % 3: Total size of the structure (Scint + Other) eqaul to Total size
                Aeq = [ones(1,2*pairs(p))];
                beq = [total_size];      
            else                    % 4: Total scint equal to Total size
                Aeq = [even];
                beq = [total_size];
            end
            
            lb = lower_bound*ones(1, 2*pairs(p));
            ub = upper_bound*ones(1, 2*pairs(p));
            
            % Randomizing starting point to improve optimal point.
            % Done several times and each time saving the best result.
            for i = 1:random_iterations
                
                d = rand(1,2 * pairs(p));
                
                % applying total size constraint
                if constraint == 1
                    d(i_scint - 1) = total_size * d(i_scint - 1) / sum(d(i_scint - 1));
                    d(i_other - 1) = total_size * d(i_other - 1) / sum(d(i_other - 1));
                elseif constraint == 2
                    d = total_size * d;
                elseif constraint == 3
                    d = total_size * d / sum(d);
                else
                    d(i_scint - 1) = total_size * d(i_scint - 1) / sum(d(i_scint - 1));
                    d(i_other - 1) = total_size * d(i_other - 1);
                end
                
                if (optimization_mode == "Efficiency")
                    % In this mode we optimize the efficiency of the
                    % structure, which is the overall emission (in all
                    % angles).
                    theta_optimization = linspace(0, pi/2, 101);
                    Scint_City_aux = @(d)-(sin(theta_optimization) * (Scint_City_fun(lambda,theta_optimization,d,n,i_scint,    coupled,control) * (Y'))) ...
                                        / (sin(theta_optimization) * (Scint_City_fun(lambda,theta_optimization,sum(d),n_bulk,2,coupled,control) * (Y'))) ;
                else % Directionality mode
                    Scint_City_aux = @(d)-sum(Scint_City_fun(lambda,theta,d,n,i_scint,coupled,control) .* (Y.'));
                end

                fprintf(['Pairs: ', num2str(pairs(p)), ', Iteration: ', num2str(i), '/', num2str(random_iterations), '\n']);
                % Calling fmincon 
                [optimized_d,F_max_val,exitflag1,output_max] = fmincon(Scint_City_aux, d, A, b, Aeq, beq, lb, ub);
                F_max_val = -F_max_val;

                if F_max_val > F_max_max
                    F_max_max = F_max_val;
                    y_max_max = optimized_d;
                    pairs_max = pairs(p);
                end
            end
        end
        
        % chosen parameters from optimization
        optimized_d = y_max_max;
        pairs = pairs_max;
        n       = [n_substrate, repmat([n_scint, n_other], 1 , pairs), 1];
        i_scint = 2:2:2 * pairs;
        
        Scint_City_aux = @(d)sum(Scint_City_fun(lambda,theta,d,n,i_scint,coupled,control) ...
                             ./ (Scint_City_fun(lambda,theta,sum(d),n_bulk,2,coupled,control)) .*(Y.'));
                    
        [y_min,F_min_val,exitflag2,output_min] = fmincon(Scint_City_aux, d, A, b, Aeq, beq, lb, ub);
    end
    
    % in known structure mode, the thicknesses are known
    if ((optimize == 0) & ~(isempty(d0)) & ~control.load_data_from_dir)
        optimized_d = d0;
    end
     
    if save_fig
        save(dir_name + "\optimized_d.mat", 'optimized_d');
        save(dir_name + "\n.mat", 'n');
        save(dir_name + "\i_scint.mat", 'i_scint');
    end
      
    % At the end of this section we calculated the thicknesses of the
    % structure. Now we will plot the efficienc and emission rate of the
    % structure.



    %% Calculting results
    %% Plotting optimal result - Efficiency
    
    % Calculating the efficiency - the emission rate outside the device,
    % devided by the emission rate inside
    
    lambda  = lambda_orig;
    theta   = theta_orig;
    if constraint == 3 % Normalizing with the Total size of the structure
        total_size = sum(optimized_d);
    else % Normalizing only with the size of the scintillator
        total_size = sum(optimized_d(i_scint - 1));
    end

if length(theta_orig) == 1
    Y      = Y_orig;
    
    Pf_max  = Scint_City_fun(lambda,theta,optimized_d,n,i_scint,inside,control);
    if ~((optimize == 0) | ~isempty(d0))
        Pf_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,inside,control);
    end
    Pf_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,inside,control);

    F_max_lambda  = Scint_City_fun(lambda,theta,optimized_d,n,i_scint,coupled,control);
    if ~((optimize == 0) | ~isempty(d0))
        F_min_lambda  = Scint_City_fun(lambda,theta,y_min,n,i_scint,coupled,control);
    end
    F_bulk_lambda = Scint_City_fun(lambda,theta,total_size,n_bulk,2,coupled,control);

    eta_max  = F_max_lambda  ./ Pf_max;
    if ~((optimize == 0) | ~isempty(d0))
        eta_min  = F_min_lambda  ./ Pf_min;
    end

    eta_bulk = F_bulk_lambda ./ Pf_bulk;
    eta_norm = max(eta_bulk .* Y);

    figure();
    plot(lambda, eta_max  .* Y  / eta_norm, 'Color', Green ,'DisplayName', '$Max$');  hold on;
    if ~((optimize == 0) | ~isempty(d0))
        plot(lambda, eta_min  .* Y  / eta_norm, 'Color', Red	 ,'DisplayName', '$Min$');  hold on;
    end
    plot(lambda, eta_bulk .* Y  / eta_norm, 'Color', Gold	 ,'DisplayName', '$bulk$'); hold on;

    graphParams(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers'], '$\lambda \ [nm]$', '$\eta (\lambda) Y(\lambda)$', 'efficiency', save_fig, dir_name); 
end
    %% Plotting optimal result - theta dependency
    
    % Plotting the emission rate as a function of the angle
    
    theta = linspace(-pi/2, pi/2, num_theta_plot);
    lambda = lambda_orig;
    F_max  = Scint_City_fun(lambda,theta,optimized_d,n,i_scint,coupled,control) * Y';
    if ~((optimize == 0) | ~isempty(d0))
        F_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,coupled,control) * Y';
    end

    F_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,coupled,control) * Y';
    norm_F = F_bulk(ceil(end / 2)); % at theta = 0

    figure();
    plot(theta, F_max  ./ norm_F,  'Color', Green, 'DisplayName', '$Max$');  hold on;
    if ~((optimize == 0) | ~isempty(d0))
        plot(theta, F_min  ./ norm_F,  'Color', Red	,  'DisplayName', '$Min$');  hold on;
    end
    plot(theta, F_bulk ./ norm_F, 'Color', Gold	,  'DisplayName', '$bulk$'); hold on;
    
    graphParams(['Optimal theoretical result - ', num2str(length(n)), ' layers'], '$\theta\ [rad]$', '$f$', 'theta', save_fig, dir_name);
    
if length(lambda_orig) > 1 & length(theta_orig) > 1
    theta = linspace(-pi/2, pi/2, num_theta_plot);
    lambda = lambda_orig;
    F_max  = Scint_City_fun(lambda,theta,optimized_d,n,i_scint,coupled,control);
    if ~((optimize == 0) | ~isempty(d0))
        F_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,coupled,control);
    end

    F_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,coupled,control);
    norm_F = F_bulk(ceil(end / 2)); % at theta = 0

    theta  = repmat(theta', 1, length(lambda));
    lambda = repmat(lambda, length(theta), 1);
    figure();
    contourf(theta, lambda, F_max  ./ norm_F,  'Color', Green, 'DisplayName', '$Max$');  hold on;
    graphParams(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers'], '$\theta\ [rad]$', '$f$', 'theta_lambda_max', save_fig, dir_name); colorbar;
    if ~((optimize == 0) | ~isempty(d0))
        figure();
        contourf(theta, lambda, F_min  ./ norm_F,  'Color', Red	,  'DisplayName', '$Min$');  hold on;
        graphParams(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers'], '$\theta\ [rad]$', '$f$', 'theta_lambda_min', save_fig, dir_name); colorbar;
    end
    figure();
    contourf(theta, lambda, F_bulk ./ norm_F, 'Color', Gold	,  'DisplayName', '$bulk$'); hold on;
    graphParams(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers'], '$\theta\ [rad]$', '$f$', 'theta_lambda_bulk', save_fig, dir_name); colorbar;
end
    %% Plotting optimized_d
    figure();
    plot(1:length(optimized_d), optimized_d, 'DisplayName', 'All');  hold on;
    plot(i_scint-1, optimized_d(i_scint - 1), 'ro', 'DisplayName', 'Scint');  hold on;
    graphParams(['Optimal thicknesses - ', num2str(2*pairs + 2), ' layers'], '$Layer\ Number$', '$Y_{max}$[nm]', 'optimized_d', save_fig, dir_name);

    %% Plotting MTF and Test image
    
    if ~ImageProcessing_Params.Show_Image_Processing
        return;
    end
    ImageProcessing_Params.DR = 1; % Dynamic range
    ImageProcessing_Params.max_distance = 1000;
    nbins = ImageProcessing_Params.nbins;
    max_distance = ImageProcessing_Params.max_distance;
    
    % Calculating psf
    [x_opt, psf_opt]   = psf_computation(optimized_d, n,      i_scint, mu_lambda,control,ImageProcessing_Params);
    [x_bulk, psf_bulk] = psf_computation(total_size,  n_bulk, 2,       mu_lambda,control,ImageProcessing_Params);
    
    image = imread('barbara.png');
    [psf_opt_2D, p_im_opt,mse_opt,MTF_opt] = Image_Processing(image,psf_opt ,sum(optimized_d),ImageProcessing_Params);
    [psf_bul_2D, p_im_bul,mse_bul,MTF_bul] = Image_Processing(image,psf_bulk,total_size,ImageProcessing_Params);
    
	[X_opt_1 , X_opt_2] = meshgrid(x_opt*1e-3, x_opt*1e-3);
    [X_bulk_1, X_bulk_2] = meshgrid(x_bulk*1e-3, x_bulk*1e-3);
    figure(); mesh(X_opt_1 , X_opt_2 , psf_opt_2D, 'DisplayName', '$Optimized$');  hold on;
    graphParams2(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers'], '$x [\mu m]$', '$y [\mu m]$', 'psf_opt_2D', save_fig, dir_name);
    figure(); mesh(X_bulk_1, X_bulk_2, psf_bul_2D, 'DisplayName', '$Bulk$');       hold on;
    graphParams2(['Bulk theoretical result - ', num2str(2*pairs + 2), ' layers'], '$x [\mu m]$', '$y [\mu m]$', 'psf_bulk_2D', save_fig, dir_name);

    psf_normalization = sum(psf_opt);
    
    max_freq = double(1/max_distance)*double(nbins/2)/nbins * 1e6;
    spatial_freq = double(1/max_distance)*double((-(nbins/2):((nbins/2) - 1)))/nbins * 1e6; % Conversion to 1/mm (instead of 1/nm)
    figure();    
    plot(spatial_freq, MTF_opt, 'Color', Green, 'DisplayName', '$Max$');  hold on;
    plot(spatial_freq, MTF_bul, 'Color', Gold , 'DisplayName', '$bulk$'); hold on;
    yline(0.03, '--', 'DisplayName', '3$\% \ line\ (resolution)$'); hold on;
    xlim([0,max_freq]);
    graphParams(['MTF - ', num2str(2*pairs + 2), ' layers'], 'LP / mm', 'MTF (\%)', 'mtf', save_fig, dir_name);
    
    
    % Test image
    figure();
    plot(x_opt*1e-3,  psf_opt / psf_normalization , 'Color', Green, 'DisplayName', '$Max$');  hold on;
    plot(x_bulk*1e-3, psf_bulk / psf_normalization, 'Color', Gold	,  'DisplayName', '$bulk$'); hold on;
    graphParams('psf', '$x [\mu m]$', '', 'psf', save_fig, dir_name);
    
    % Plotting
    figure();
    subplot(1,3, 1); imshow(image);    graphParams2('Original Image', '', '', '', 0, '');
    subplot(1,3, 2); imshow(p_im_bul); graphParams2(['Bulk: MSE = ', num2str(mse_bul)],      '', '', '', 0, '');
    subplot(1,3, 3); imshow(p_im_opt); graphParams2(['Optimized: MSE = ', num2str(mse_opt)], '', '', 'barbara', save_fig, dir_name);
    
    
    %% Save
    if save_fig
        save([dir_name,'\workspace']);
    end
    
    %% Lumerical Simulation
    Simulate = simulation_params.simulate;
    if Simulate & ~is_Gz
        lumerical_path = simulation_params.lumerical_path;
        
        
        structure_name = 'Sim';
        span		= max(lambda_orig) - min(lambda_orig);
        num_lambda	= length(lambda_orig);

        Lumerical.file_name = ['Lumerical_',structure_name]; 
        Lumerical.mesh_accu = 3; % Affect results accuracy - Recommended value for fast simulation: 3. Recommended value for high accuracy: 7.
        Lumerical.sim_width_param = 10; % Affect high angles result accuracy. Recommended minimal value: 10.
        Lumerical.lambda = mu_lambda;
        Lumerical.span = span;
        Lumerical.num_lambda = num_lambda;
        Lumerical.d = optimized_d;
        Lumerical.dz = dz_Nz;
        Lumerical.num_theta = num_theta_plot;
        Lumerical.i_scint = i_scint;
        Lumerical.n = n;
        save('Lumerical.mat','Lumerical'); 
        command = ['"', lumerical_path, '" -run Scint_City_Sim.lsf -exit'];
%         command = ['"', lumerical_path, '" -run Scint_City_Sim.lsf'];
        return_val = system(command);
    %     assert(logical(return_val) ,'Simulation failed to run');
    
        load(['Lumerical_', structure_name], 'Lumerical_results');
        
        Lumerical_results_Opt = Lumerical_results;

        F_Sim_Opt = Lumerical_results_Opt.X+Lumerical_results_Opt.Y+Lumerical_results_Opt.Z;
        
        Lumerical.d = total_size;
        Lumerical.i_scint = 2;
        Lumerical.n = n_bulk;
        save('Lumerical.mat','Lumerical'); 
        
        return_val = system(command);
        
        load(['Lumerical_', structure_name], 'Lumerical_results');
        
        Lumerical_results_Bulk = Lumerical_results;
        
        F_Sim_Bulk = Lumerical_results_Bulk.X+Lumerical_results_Bulk.Y+Lumerical_results_Bulk.Z;
        
        c = 299792458;
        sim_lambda = c./Lumerical_results.f0' * 1e9;
        sim_theta = Lumerical_results.theta0';
        sim_d = Lumerical_results.d * 1e9;
        
        lambda = lambda_orig;
        F_max  = Scint_City_fun(lambda,sim_theta,optimized_d,n,i_scint,coupled,control);

        F_bulk = Scint_City_fun(lambda,sim_theta,total_size,n_bulk,2,coupled,control);
        norm_F_theo = sum(F_bulk(ceil(end / 2),:).*Y); % at theta = 0
        norm_F_sim = sum(F_Sim_Bulk(ceil(end / 2),:).*Y); % at theta = 0
        
        figure;
        hold on;
        plot(sim_theta, sum(F_max.*Y,2)/norm_F_theo,  'Color', Green, 'DisplayName', 'Opt - Theory');
        plot(sim_theta, sum(F_Sim_Opt.*Y,2)/norm_F_sim, ':',  'Color', Green_Sim, 'DisplayName', 'Opt - Simulation');
        plot(sim_theta, sum(F_bulk.*Y,2)/norm_F_theo,  'Color', Gold, 'DisplayName', 'Bulk - Theory');
        plot(sim_theta, sum(F_Sim_Bulk.*Y,2)/norm_F_sim, ':',  'Color', Gold_Sim, 'DisplayName', 'Bulk - Simulation');

        legend('show');
        grid on;
        xlabel('$\theta$ [rad]');
        ylabel('Angular power');
        hold off;

        % xticklabels({});
        % yticklabels({});
        xticks([-pi/2 -pi/4 0 pi/4 pi/2]);
        % yticks([0 1 2 3 4 5 6]);
        set(gca,'linewidth',2.5);
        set(gca, 'FontSize', 14);
        box on;
    end
    
    %% Save
    if save_fig
        save([dir_name,'\workspace']);
    end
    
end

function graphParams(ptitle, pxlabel, pylabel, figname, save_fig, dir_name) 
    grid on;
    title(ptitle);
    xlabel(pxlabel);
    ylabel(pylabel);
    set(gca, 'FontSize', 14);
    set(gcf,'color','w');
    set(gca,'linewidth',2.5);
    legend('show');
    if save_fig
        saveas(gcf, dir_name + "\" + figname + ".svg");
        saveas(gcf, dir_name + "\" + figname + ".fig");
    end
    %set(gca,'XTickLabel',[], 'YTickLabel',[]);
end
function graphParams2(ptitle, pxlabel, pylabel, figname, save_fig, dir_name) 
    grid on;
    title(ptitle);
    xlabel(pxlabel);
    ylabel(pylabel);
    set(gca, 'FontSize', 14);
    set(gcf,'color','w');
    set(gca,'linewidth',2.5);
    if save_fig
        saveas(gcf, dir_name + "\" + figname + ".svg");
        saveas(gcf, dir_name + "\" + figname + ".fig");
    end
    %set(gca,'XTickLabel',[], 'YTickLabel',[]);
end