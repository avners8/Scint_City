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

    if ~((optimize == 0) & ~isempty(d0))
        n       = [n_si, repmat([n_scint, n_ox], 1 , pairs), 1];
        i_scint = 2:2:2 * pairs;
    end

    %% Some stuff

    theta_orig = theta;
    lambda_orig = lambda;

    Green = [0.4660 0.6740 0.1880];
    Red   = [0.6350 0.0780 0.1840];
    Gold  = [0.9290 0.6940 0.1250];

    dir_name = [pwd, '\', dir_name];

    %% Calculating constraints

    % starting point
    d =  total_size / pairs * ones(1, 2 * pairs);
    even = 0.5 * (1+(-1).^(1:length(d)));
    odd  = 0.5 * (1+(-1).^(0:length(d) - 1));

    % Defining the constraints - Ad<=b
    % Total oxide size is equal to 1 micron 
    A = [];
    b = [];
    Aeq = [even;odd];
    beq = [total_size,total_size];
    lb = zeros(1, length(d));
    ub = total_size*ones(1, length(d));
    % Total oxide size can be smaller than 1 micron 
    % A = [ -odd; odd; -diag(ones(1, length(d)))];
    % b = [ -micron, micron, -dz*ones(1, length(d))];


    %% Calculating the minimum using Matlab Optimization Toolbox

    coupled = true;
    inside  = false;
	   
    if sigma_lambda == 0
       lambda = mu_lambda;
       Y      = 1;
    end
    n_bulk  = [n_si n_scint 1];
    i_other = zeros(1,length(n) - length(i_scint) - 2);
    count_scint = 1; count_other = 1;
    for i = 2:length(n) - 1
        if i - 1 == i_scint(count_scint)
            count_scint = count_scint + 1;
        else
            i_other(count_other) = i - 1;
            count_other = count_other + 1;
        end
    end
    
    if optimize
        y_max_max = d;
        F_max_max = 0;
        for i = 1:random_iterations
            d = rand(size(d));
            
            if total_size ~= 0
                d(i_scint) = total_size * d(i_scint) / sum(d(i_scint));
                d(i_other) = total_size * d(i_other) / sum(d(i_other));
            end
            
            Scint_City_aux = @(d)-sum(Scint_City_fun(lambda,theta,d,n,i_scint,dz_Nz,coupled,is_Gz) ...
                                  ./ (Scint_City_fun(lambda,theta,sum(d),n_bulk,2,dz_Nz,coupled,is_Gz)) .*(Y.'));
                              
            [y_max,F_max_val,exitflag1,output_max] = fmincon(Scint_City_aux, d, A, b, Aeq, beq, lb, ub);
            F_max_val = -F_max_val;
            
            if F_max_val > F_max_max
                F_max_max = F_max_val;
                y_max_max = y_max;
            end
        end
        y_max = y_max_max;
        
        Scint_City_aux = @(d)sum(Scint_City_fun(lambda,theta,d,n,i_scint,dz_Nz,coupled,is_Gz) ...
                             ./ (Scint_City_fun(lambda,theta,sum(d),n_bulk,2,dz_Nz,coupled,is_Gz)) .*(Y.'));
                         
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
    total_size = sum(y_max(i_scint));
    
    Pf_max  = Scint_City_fun(lambda,theta,y_max,n,i_scint,dz_Nz,inside,is_Gz);
    if ~((optimize == 0) & ~isempty(d0))
        Pf_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,dz_Nz,inside,is_Gz);
    end
    Pf_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,dz_Nz,inside,is_Gz);

    F_max_lambda  = Scint_City_fun(lambda,theta,y_max,n,i_scint,dz_Nz,coupled,is_Gz);
    if ~((optimize == 0) & ~isempty(d0))
        F_min_lambda  = Scint_City_fun(lambda,theta,y_min,n,i_scint,dz_Nz,coupled,is_Gz);
    end
    F_bulk_lambda = Scint_City_fun(lambda,theta,total_size,n_bulk,2,dz_Nz,coupled,is_Gz);

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

    theta = linspace(-pi/2, pi/2, 301);
    lambda = central_lambda;
    F_max  = Scint_City_fun(lambda,theta,y_max,n,i_scint,dz_Nz,coupled,is_Gz);
    if ~((optimize == 0) & ~isempty(d0))
        F_min  = Scint_City_fun(lambda,theta,y_min,n,i_scint,dz_Nz,coupled,is_Gz);
    end

    F_bulk = Scint_City_fun(lambda,theta,total_size,n_bulk,2,dz_Nz,coupled,is_Gz);
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
    plot(1:length(y_max), y_max * 1e9);  hold on;
    graphParams(['Optimal thicknesses - ', num2str(2*pairs + 2), ' layers, ', '$lambda=$', num2str(lambda*1e9), '[nm]'], '$Layer\ Number$', '$Y_{max}$[nm]', 'y_max', save_fig, dir_name);
    y_max * 1e9


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
