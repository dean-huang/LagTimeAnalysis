% Find the slopes for an n leg piecewise linear plot, given arbitrary 
% x values for control points.

% Guess structure: Slot 1 is the y-intercept and slots 2 to (n+1) are the 
% y value slope guesses. Slots (n+2) to (2n+2) are the x positions of the 
% control points, including endpoints. Note that if n is the number of 
% legs, n+1 would be the number of control points.

function [fitparam,param_err,yfit,Jacob] = FitPiecewiseLinear_nLeg_Slopes( x, y_exp, guess, debug, sigma )

    if ~exist('debug','var')
        debug = false;
    end
    if ~exist('sigma','var')
        sigma = sqrt(mean((y_exp(2:end)-y_exp(1:end-1)).^2,'omitnan'))/sqrt(2);
    end

    N = numel(x);
    
    n_ctrl = numel(guess)/2;
    param0 = guess(1:n_ctrl);
    xpts = guess(n_ctrl+1:end);
    
    [fitparam,~,~,~,~,~,Jacob] = lsqnonlin( @fitdiff, param0 );
    yfit = Theory(x,xpts,fitparam,N);
    
    Info = full((Jacob'*Jacob)/(sigma^2));
    param_err = sqrt(diag(inv(Info)))';
    
    function dy = fitdiff( param )
        
        y_th = Theory(x,xpts,param,N);
        dy = ( y_exp - y_th );
        
        if debug
            
            figure(1);
            clf;
        
            plot( x, y_exp, 'r.' );
            hold on;
            
            plot( x, y_th, 'b.-' );
            
            drawnow;
            pause(.01);
            
        end
    end
end

function y = Theory(x,xpts,param,N)
    % xpts is horizontal value of control points.
    
    y = zeros(size(x));
    ypts = zeros(size(param));
    ypts(1) = param(1); % Recall that the first element is the y-intercept.
    for ii = 2:numel(ypts)
        ypts(ii) = ypts(ii-1)+param(ii)*(xpts(ii)-xpts(ii-1));
    end
    % ypts now represents y values of the control points.
    
    for ii = 1:N
        xx = x(ii);
        jj = sum(xpts<=xx);
        if jj > 0
            if jj < numel(ypts)
                xstart = xpts(jj);
                xend = xpts(jj+1);
                ystart = ypts(jj);
                yend = ypts(jj+1);
                y(ii) = ystart + ((yend-ystart)/(xend-xstart))*(xx-xstart);
            else
                y(ii) = ypts(jj);
            end
        end
    end
    

    
    
end
