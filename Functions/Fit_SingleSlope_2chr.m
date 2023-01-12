% Given data from 2chr Vibrio, fit it all with a single slope value. The
% legs will be connected for each chromosome, but not between the two chrs.

% Guess structure: Slots 1 to (n+1) are the x control points, including 
% endpoints.

function [fitparam,param_err,yfit] = Fit_SingleSlope_2chr( x, y_exp, xpts, sigma, debug )

    if ~exist('sigma','var') || isempty(sigma)
        sigma = sqrt(mean((y_exp(2:end)-y_exp(1:end-1)).^2,'omitnan'))/sqrt(2);
    end
    if ~exist('debug','var') || isempty(debug)
        debug = false;
    end

    param0 = zeros(3,1);
    
    keeppts = x >= min(xpts) & x <= max(xpts);
    x = x(keeppts);
    y_exp = y_exp(keeppts);
    
    N = numel(x);
    
    % Display off turns off the optimization complete statement. Lower and
    % upper bounds are empty matrices.
    opts = optimset('Display','off');
    [fitparam,~,~,~,~,~,Jacob] = lsqnonlin( @fitdiff, param0, [],[], opts);
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
    ystarts = param(1:2);
    m = param(end);
    
    for ii = 1:N
        xx = x(ii);
        jj = sum(xpts<=xx);
        xstart = xpts(jj);
        if jj < numel(xpts)
            if mod(jj,2) == 1
                ystart = ystarts(ceil(jj/2));
            else
                ystart = ystarts(jj/2)+m*(xpts(jj)-xpts(jj-1));
            end
            updown = (-1).^(mod(jj,2)+1);
            y(ii) = ystart + updown*m*(xx-xstart);
        else
            y(ii) = ystarts(2);
        end
    end
    
end
