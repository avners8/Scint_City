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
                                total_size_con, ...
                                optimize,       ...
                                save_fig,       ...
                                dir_name )
    
    set(groot, 'defaultFigurePosition', [100 100 900 600]); % figure size
    set(groot, 'defaultTextInterpreter', 'latex'); % latex
    set(groot, 'defaultLegendInterpreter', 'latex'); % latex
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); % latex
    set(0, 'DefaultLineLineWidth', 4);

    %% Building the structure
    
    if sigma_lambda == 0
        Y_orig = zeros(size(lambda));
        Y_orig(lambda == mu_lambda) = 1;
    else
        Y_orig = normpdf(lambda, mu_lambda, sigma_lambda);
    end
    Y      = Y_orig;
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
    n_bulk  = [n_si n_scint 1];
    
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
            if total_size_con
                Aeq = [even];
                beq = [total_size];
%                 Aeq = [even;odd];
%                 beq = [total_size,total_size];
            else
                Aeq = [];
                beq = [];
            end
            lb = zeros(1, 2*pairs(p));
            ub = total_size*ones(1, 2*pairs(p));
            
            % Randomizing starting point to improve optimal point
            for i = 1:random_iterations
                d = rand(1,2 * pairs(p));
                
                % applying total size constraint
                if total_size_con
                    d(i_scint) = total_size * d(i_scint) / sum(d(i_scint));
                    d(i_other) = total_size * d(i_other) / sum(d(i_other));
                else
                    d = total_size * d;
                end

                Scint_City_aux = @(d)-sum(Scint_City_fun(lambda,theta,d,n,i_scint,coupled,control) ...
                                      ./ (Scint_City_fun(lambda,theta,sum(d),n_bulk,2,coupled,control)) .*(Y.'));

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
        
        d = rand(1,2 * pairs);            
        [y_min,F_min_val,exitflag2,output_min] = fmincon(Scint_City_aux, d, A, b, Aeq, beq, lb, ub);
    end
    
    if ((optimize == 0) & ~(isempty(d0)))
        y_max = d0;
    end
    
    %% Calculting results
    %% Plotting optimal result - Efficiency

    Y       = Y_orig;
    lambda  = lambda_orig;
    theta   = theta_orig;
    total_size = sum(y_max(i_scint - 1));
    
    Pf_max  = Scint_City_fun(lambda,theta,y_max,n,i_scint,inside,control);
    if ~((optimize == 0) & ~isempty(d0))
        Pf_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,inside,control);
    end
    Pf_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,inside,control);

    F_max_lambda  = Scint_City_fun(lambda,theta,y_max,n,i_scint,coupled,control);
    if ~((optimize == 0) & ~isempty(d0))
        F_min_lambda  = Scint_City_fun(lambda,theta,y_min,n,i_scint,coupled,control);
    end
    F_bulk_lambda = Scint_City_fun(lambda,theta,total_size,n_bulk,2,coupled,control);

    eta_max  = F_max_lambda  ./ Pf_max;
    if ~((optimize == 0) & ~isempty(d0))
        eta_min  = F_min_lambda  ./ Pf_min;
    end

    eta_bulk = F_bulk_lambda ./ Pf_bulk;
    eta_norm = max(eta_bulk.' .* Y);

    figure();
    plot(lambda*1e9, eta_max.'  .* Y  / eta_norm, 'Color', Green ,'DisplayName', '$Max$');  hold on;
    if ~((optimize == 0) & ~isempty(d0))
        plot(lambda*1e9, eta_min.'  .* Y  / eta_norm, 'Color', Red	 ,'DisplayName', '$Min$');  hold on;
    end
    plot(lambda*1e9, eta_bulk.' .* Y  / eta_norm, 'Color', Gold	 ,'DisplayName', '$bulk$'); hold on;

    graphParams(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers'], '$\lambda \ [nm]$', '$\eta (\lambda) Y(\lambda)$', 'efficiency', save_fig, dir_name); 

    %% Plotting optimal result - theta dependency

    theta = linspace(-pi/2, pi/2, 1001);
    lambda = central_lambda;
    F_max  = Scint_City_fun(lambda,theta,y_max,n,i_scint,coupled,control);
    if ~((optimize == 0) & ~isempty(d0))
        F_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,coupled,control);
    end

    F_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,coupled,control);
    norm_F = F_bulk(ceil(end / 2)); % at theta = 0

    figure();
    plot(theta, F_max  ./ norm_F,  'Color', Green, 'DisplayName', '$Max$');  hold on;
    if ~((optimize == 0) & ~isempty(d0))
        plot(theta, F_min  ./ norm_F,  'Color', Red	,  'DisplayName', '$Min$');  hold on;
    end
    plot(theta, F_bulk ./ norm_F, 'Color', Gold	,  'DisplayName', '$bulk$'); hold on;
    
    graphParams(['Optimal theoretical result - ', num2str(2*pairs + 2), ' layers, ', '$lambda=$', num2str(lambda*1e9), '[nm]'], '$\theta\ [rad]$', '$f$', 'theta', save_fig, dir_name);
       
    %% Plotting y_max
    figure();
    plot(1:length(y_max), y_max * 1e9, 'DisplayName', 'All');  hold on;
    plot(i_scint-1, y_max(i_scint - 1) * 1e9, 'ro', 'DisplayName', 'Scint');  hold on;
    graphParams(['Optimal thicknesses - ', num2str(2*pairs + 2), ' layers, ', '$lambda=$', num2str(lambda*1e9), '[nm]'], '$Layer\ Number$', '$Y_{max}$[nm]', 'y_max', save_fig, dir_name);
    y_max_nm = y_max * 1e9;
    if save_fig
        save(dir_name + "\y_max.mat", 'y_max_nm');
    end
%     save('workspace.mat', '-append');
    %% Plotting MTF and Test image
%     load('workspace.mat');
    
    ImageProcessing_Params.nbins = 256;
    ImageProcessing_Params.h     = 0;
    ImageProcessing_Params.DR = 1; % Dynamic range
    ImageProcessing_Params.sigma_onion_noise = 1;
    
    % Calculating psf
    [x_opt, psf_opt]   = psf_computation(y_max,      n,      i_scint, central_lambda,control,ImageProcessing_Params);
    [x_bulk, psf_bulk] = psf_computation(total_size, n_bulk, 2,       central_lambda,control,ImageProcessing_Params);
    
    image = imread('barbara.png');
    [p_im_opt,mse_opt,MTF_opt] = Image_Processing(image,psf_opt ,sum(y_max),ImageProcessing_Params);
    [p_im_bul,mse_bul,MTF_bul] = Image_Processing(image,psf_bulk,total_size,ImageProcessing_Params);
    
    figure();    
    plot(ceil(length(MTF_opt)/2):length(MTF_opt), MTF_opt(ceil(length(MTF_opt)/2) :length(MTF_opt)), 'Color', Green, 'DisplayName', '$Max$');  hold on;
    plot(ceil(length(MTF_bul)/2):length(MTF_bul), MTF_bul(ceil(length(MTF_bul)/2) :length(MTF_bul)), 'Color', Gold , 'DisplayName', '$bulk$'); hold on;
    xlim([ceil(length(MTF_bul)/2)+1,length(MTF_bul)]);
    graphParams(['MTF - ', num2str(2*pairs + 2), ' layers'], '$Spacial\ Frequency$', '$log(1+MTF)$', 'mtf', save_fig, dir_name);
    
    
    % Test image
    figure();
    plot(x_opt,  psf_opt , 'Color', Green, 'DisplayName', '$Max$');  hold on;
    plot(x_bulk, psf_bulk, 'Color', Gold	,  'DisplayName', '$bulk$'); hold on;
    graphParams('psf', 'pixel', '', 'psf', save_fig, dir_name);
    
    % Plotting
    figure();
    subplot(1,3, 1); imshow(image);    graphParams2('Original Image', '', '', '', 0, '');
    subplot(1,3, 2); imshow(p_im_bul); graphParams2(['Bulk: MSE = ', num2str(mse_bul)],      '', '', '', 0, '');
    subplot(1,3, 3); imshow(p_im_opt); graphParams2(['Optimized: MSE = ', num2str(mse_opt)], '', '', 'barbara', save_fig, dir_name);
    
    
end

function graphParams(ptitle, pxlabel, pylabel, figname, save_fig, dir_name) 
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
    legend('show');
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