function One_D_PhC_Optimization(dz_Nz,          ...
                                theta,          ...
                                lambda,         ...
                                central_lambda, ...
                                mu_lambda,      ...
                                sigma_lambda,   ...
                                pairs,          ...
                                n_si,           ...
                                n_scint,        ...
                                n_ox,           ...
                                d0,             ...
                                n,              ...
                                i_scint,        ...
                                is_Gz,          ...
                                random_iterations, ...
                                total_size,     ...
                                constraint, ...
                                optimize,       ...
                                ImageProcessing_Params, ...
                                save_fig,       ...
                                dir_name )
    %% Building the structure
    
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
    control.dz_Nz = dz_Nz; control.is_Gz = is_Gz;  control.sum_on_z = true;

    %% Some stuff

    theta_orig = theta;
    lambda_orig = lambda;

    Green = [0.4660 0.6740 0.1880];
    Red   = [0.6350 0.0780 0.1840];
    Gold  = [0.9290 0.6940 0.1250];

    dir_name = [pwd, '\', dir_name];


    %% Calculating the minimum using Matlab Optimization Toolbox

    coupled = true;
    inside  = false;
	   
    if sigma_lambda == 0
       lambda = mu_lambda;
       Y      = 1;
    end
    if optimize
        n_bulk  = [n_si n_scint 1];
    else
        n_bulk  = [n(1) n(i_scint(1)) n(end)];
    end
    
    if optimize
        y_max_max = 0;
        F_max_max = 0;
        pairs_max = pairs(1);
        
        % Iterating on number of layers in the structure
        for p = 1:length(pairs)
            
            n       = [n_si, repmat([n_scint, n_ox], 1 , pairs(p)), 1];
            i_scint = 2:2:2 * pairs(p);
            
            i_other = zeros(1,length(n) - length(i_scint) - 2);
            count_scint = 1; count_other = 1;
            for j = 2:length(n) - 1
                if j - 1 == i_scint(count_scint)
                    count_scint = count_scint + 1;
                else
                    i_other(count_other) = j - 1;
                    count_other = count_other + 1;
                end
            end
            
            % Calculating constraints
            % starting point
            even = 0.5 * (1+(-1).^(1:(2*pairs(p))));
            odd  = 0.5 * (1+(-1).^(0:(2*pairs(p) - 1)));

            % Defining the constraints - Ad=b
            % Total oxide size is equal to total size 
            A = [];
            b = [];
            if constraint == 1      % 1: Total scint and Total other equal Total size
                Aeq = [even];
                beq = [total_size];
%                 Aeq = [even;odd];
%                 beq = [total_size,total_size];
            elseif constraint == 2  % 2: Each layer is between 0 and Total size
                Aeq = [];
                beq = [];
            else                    % 3: Total size of the structure (Scint + Other) eqaul to Total size
                Aeq = [ones(1,2*pairs(p))];
                beq = [total_size];                
            end
            lb = zeros(1, 2*pairs(p));
            ub = total_size*ones(1, 2*pairs(p));
            
            % Randomizing starting point to improve optimal point
            for i = 1:random_iterations
                
                d = rand(1,2 * pairs(p));
                
                % applying total size constraint
                if constraint == 1
                    d(i_scint) = total_size * d(i_scint) / sum(d(i_scint));
                    d(i_other) = total_size * d(i_other) / sum(d(i_other));
                elseif constraint == 2
                    d = total_size * d;
                else
                    d = total_size * d / sum(d);
                end

                Scint_City_aux = @(d)-sum(Scint_City_fun(lambda,theta,d,n,i_scint,coupled,control) ...
                                      ./ (Scint_City_fun(lambda,theta,sum(d),n_bulk,2,coupled,control)) .*(Y.'));

                fprintf(['Pairs: ', num2str(pairs(p)), ', Iteration: ', num2str(i), '/', num2str(random_iterations), '\n']);
                % Calling fmincon 
                [y_max,F_max_val,exitflag1,output_max] = fmincon(Scint_City_aux, d, A, b, Aeq, beq, lb, ub);
                F_max_val = -F_max_val;

                if F_max_val > F_max_max
                    F_max_max = F_max_val;
                    y_max_max = y_max;
                    pairs_max = pairs(p);
                end
            end
        end
        
        % chosen parameters from optimization
        y_max = y_max_max;
        pairs = pairs_max;
        n       = [n_si, repmat([n_scint, n_ox], 1 , pairs), 1];
        i_scint = 2:2:2 * pairs;
        
        Scint_City_aux = @(d)sum(Scint_City_fun(lambda,theta,d,n,i_scint,coupled,control) ...
                             ./ (Scint_City_fun(lambda,theta,sum(d),n_bulk,2,coupled,control)) .*(Y.'));
                    
        [y_min,F_min_val,exitflag2,output_min] = fmincon(Scint_City_aux, d, A, b, Aeq, beq, lb, ub);
    end
    
    if ((optimize == 0) & ~(isempty(d0)))
        y_max = d0;
    end
    
    %% Calculting results
    %% Plotting optimal result - Efficiency
    lambda  = lambda_orig;
    theta   = theta_orig;
    if constraint == 3 % Normalizing with the Total size of the structure
        total_size = sum(y_max);
    else % Normalizing only with the size of the scintillator
        total_size = sum(y_max(i_scint - 1));
    end

if length(theta_orig) == 1
    Y      = Y_orig;
    
    Pf_max  = Scint_City_fun(lambda,theta,y_max,n,i_scint,inside,control);
    if ~((optimize == 0) | ~isempty(d0))
        Pf_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,inside,control);
    end
    Pf_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,inside,control);

    F_max_lambda  = Scint_City_fun(lambda,theta,y_max,n,i_scint,coupled,control);
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

    theta = linspace(-pi/2, pi/2, 1001);
    lambda = central_lambda;
    F_max  = Scint_City_fun(lambda,theta,y_max,n,i_scint,coupled,control);
    if ~((optimize == 0) | ~isempty(d0))
        F_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,coupled,control);
    end

    F_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,coupled,control);
    norm_F = F_bulk(ceil(end / 2)); % at theta = 0

    figure();
    plot(theta, F_max  ./ norm_F,  'Color', Green, 'DisplayName', '$Max$');  hold on;
    if ~((optimize == 0) | ~isempty(d0))
        plot(theta, F_min  ./ norm_F,  'Color', Red	,  'DisplayName', '$Min$');  hold on;
    end
    plot(theta, F_bulk ./ norm_F, 'Color', Gold	,  'DisplayName', '$bulk$'); hold on;
    
    graphParams(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers, ', '$lambda=$', num2str(lambda), '[nm]'], '$\theta\ [rad]$', '$f$', 'theta', save_fig, dir_name);
    
if length(lambda_orig) > 1 & length(theta_orig) > 1
    theta = linspace(-pi/2, pi/2, 1001);
    lambda = lambda_orig;
    F_max  = Scint_City_fun(lambda,theta,y_max,n,i_scint,coupled,control);
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
    %% Plotting y_max
    figure();
    plot(1:length(y_max), y_max, 'DisplayName', 'All');  hold on;
    plot(i_scint-1, y_max(i_scint - 1), 'ro', 'DisplayName', 'Scint');  hold on;
    graphParams(['Optimal thicknesses - ', num2str(2*pairs + 2), ' layers'], '$Layer\ Number$', '$Y_{max}$[nm]', 'y_max', save_fig, dir_name);
    y_max_nm = y_max;
    if save_fig
        save(dir_name + "\y_max.mat", 'y_max_nm');
    end

    %% Plotting MTF and Test image
    
    ImageProcessing_Params.DR = 1; % Dynamic range
    ImageProcessing_Params.max_distance = 1000;
    nbins = ImageProcessing_Params.nbins;
    max_distance = ImageProcessing_Params.max_distance;
    
    % Calculating psf
    [x_opt, psf_opt]   = psf_computation(y_max,      n,      i_scint, central_lambda,control,ImageProcessing_Params);
    [x_bulk, psf_bulk] = psf_computation(total_size, n_bulk, 2,       central_lambda,control,ImageProcessing_Params);
    
    image = imread('barbara.png');
    [p_im_opt,mse_opt,MTF_opt] = Image_Processing(image,psf_opt ,sum(y_max),ImageProcessing_Params);
    [p_im_bul,mse_bul,MTF_bul] = Image_Processing(image,psf_bulk,total_size,ImageProcessing_Params);
    
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