% Fit a sigmoid function to data. 

% f(x) = c + L/(1+exp(-k(x-x_0)))

% Parameters:
% (1) Growth rate, k. (2) Midpoint, x_0. (3) Capacity, L. (4) Lower limit, c. 

function [fitparam,param_err,yfit] = Fit_Sigmoid( x, y_expt, guess, lb, ub, sigma, debug )

    if ~exist('debug','var') || isempty(debug)
        debug = false;
    end

    param0 = guess;
    
    % Display off turns off the optimization complete statement. Lower and
    % upper bounds are empty matrices.
    opts = optimset('Display','off');
    [fitparam,~,~,~,~,~,Jacob] = lsqnonlin( @fitdiff, param0, lb, ub, opts);
    yfit = Theory(x,fitparam);
    
    if ~exist('sigma','var') || isempty(sigma)
        sigma = sqrt(mean(fitdiff(fitparam).^2));
    end

    Info = full((Jacob'*Jacob)/(sigma^2));
    param_err = sqrt(diag(inv(Info)))';
    
    function dy = fitdiff( param )
        
        y_th = Theory(x,param);
        dy = ( y_expt - y_th );
        
        if debug
            
            figure(1);
            clf;
        
            plot( x, y_expt, 'r.', 'MarkerSize', 20 );
            hold on;
            
            xfit = (1:1000)/1000*(max(x)-min(x))+min(x);
            y_thfit = Theory(xfit,param);
            plot( xfit, y_thfit, 'b.-', 'MarkerSize', 1 );

            plot( x, y_th, 'g.-' );
            
            drawnow;
            pause(.01);
            
        end
    end
end

function y = Theory(x,param)

    k = param(1);
    x_0 = param(2);
    L = param(3);
    c = param(4);
    
    y = c + L./(1+exp(-k.*(x-x_0)));
    
end