% Find the y values of the control points for an n leg piecewise linear
% plot, given arbitrary x values for control points.

% Guess structure: Slots 1 to (n+1) are the y value guesses, including 
% endpoints. Slots (n+2) to (2n+2) are the x positions of the control 
% points, including endpoints. Note that if n is the number of legs,
% n+1 would be the number of control points.

function [fitparam,param_err,yfit,rep,rep_err] = Fit_PiecewiseLinear_nLeg_SetX( x, y_exp, guess, kG, sigma, debug )

    if ~exist('kG','var') || isempty(kG)
        disp('NEED TO INPUT kG')
        kG = 1000;
    end
    if ~exist('sigma','var') || isempty(sigma)
        sigma = Get_SigmaEst(y_exp);
    end
    if ~exist('debug','var') || isempty(debug)
        debug = false;
    end

    n_ctrl = numel(guess)/2;
    param0 = guess(1:n_ctrl);
    xpts = guess(n_ctrl+1:end);
    
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
    
    [divdelL,ydif,Slopes2Vels] = deal(zeros(numel(xpts)));
    divdelL(1,1) = 1;
    ydif(1,1) = 1;
    for ii = 2:numel(xpts)
        divdelL(ii,ii) = 1/(xpts(ii)-xpts(ii-1));
        ydif(ii,ii) = 1;
        ydif(ii,ii-1) = -1;
    end
    CtrlPts2Slopes = divdelL*ydif;
    slope = CtrlPts2Slopes*fitparam';
    slopeInfoinv = CtrlPts2Slopes*(Info\CtrlPts2Slopes');
    %slope_err = sqrt(diag(slopeInfoinv))';
    %slope = slope(2:end)';
    %slope_err = slope_err(2:end);
    
    Slopes2Vels(1,1) = 1;
    for ii = 2:numel(xpts)
        Slopes2Vels(ii,ii) = -kG.*(slope(ii).^(-2));
    end
    rep = Slopes2Vels*slope;
    % Note: slopeInfoinv is already the inverse Info matrix sandwiched by
    % two matrices, that's why here we just do regular matrix
    % multiplication instead of the left divide.
    repInfoinv = Slopes2Vels*(slopeInfoinv*Slopes2Vels');
    rep_err = sqrt(diag(repInfoinv))';
    rep = rep(2:end)';
    rep_err = rep_err(2:end);
    
    function dy = fitdiff( param )
        
        y_th = Theory(x,xpts,param,N);
        dy = ( y_exp - y_th );
        
        if debug
            
            figure(1);
            clf;
        
            plot( x, y_exp, 'r.' );
            hold on;
            
            plot( x, y_th, 'b.-' );
            
            plot( xpts, param, 'o', 'MarkerSize', 10, 'LineWidth', 3 );
            
            drawnow;
            pause(.01);
            
        end
    end
end

function y = Theory(x,xpts,param,N)
    % xpts is horizontal value of control points.
    y = zeros(size(x));
    ypts = param; 
    
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
            elseif jj == numel(ypts)
                y(ii) = ypts(end);
            else
                % Sets points that fall outside of the xpts domain to 0.
                % This makes sure that they don't affect the positioning of
                % the control points during optimization. This is for the
                % RHS.
                y(ii) = 0;
            end
        else
            % Sets points that fall outside of the xpts domain to 0 for
            % LHS.
            y(ii) = 0;
        end
    end
    
end
