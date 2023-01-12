% Fit a sine function to data. Centered around 0.

% Parameters:
% (1) Amplitude. (2) Period. (3) Phase. 

function [fitparam,param_err,yfit] = Fit_Sine( tau, y_expt, guess, lb, ub, sigma, debug )
    
    if ~exist('sigma','var') || isempty(sigma)
        sigma = Get_SigmaEst(y_expt);
    end
    if ~exist('debug','var') || isempty(debug)
        debug = false;
    end

    param0 = guess;
    
    % Display off turns off the optimization complete statement. Lower and
    % upper bounds are empty matrices.
    opts = optimset('Display','off');
    [fitparam,~,~,~,~,~,Jacob] = lsqnonlin( @fitdiff, param0, lb, ub, opts);
    yfit = Theory(tau,fitparam);
    
    Info = full((Jacob'*Jacob)/(sigma^2));
    param_err = sqrt(diag(inv(Info)))';
    
    function dy = fitdiff( param )
        
        y_th = Theory(tau,param);
        dy = ( y_expt - y_th );
        
        if debug
            
            figure(1);
            clf;
        
            plot( tau, y_expt, 'r.', 'MarkerSize', 20 );
            hold on;
            
            taufit = (1:1000)/50;
            y_thfit = Theory(taufit,param);
            plot( taufit, y_thfit, 'b.-', 'MarkerSize', 1 );

            plot( tau, y_th, 'g.-' );
            
            drawnow;
            pause(.01);
            
        end
    end
end

function y = Theory(tau,param)

    amp = param(1);
    period = param(2);
    phase = param(3);
    
    y = amp*sin(2*pi*tau/period - phase);
    
end