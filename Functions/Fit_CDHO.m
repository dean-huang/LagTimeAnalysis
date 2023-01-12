% Fit two different slopes for CD+NT (one slope) and HO regions. 

% x_bp is the x values at bp resolution. We need these to properly
% interpret the gene_CDHO positions.

% The gene_CDHO input is an array denoting which bp locations are CD+NT and
% which are HO, +1 for CD, 0 for NT, and -1 for HO. It already takes
% into account the direction of travel for the replisome, so there is a
% parity term.

% Guess structure: The first parameter is the CD+NT slope and the second is
% the HO slope. The third parameter is the starting height.

function [fitparam,param_err,yfit] = Fit_CDHO( x_bp, y_expt, gene_CDHO, guess, sigma, debug )

    if ~exist('sigma','var') || isempty(sigma)
        sigma = Get_SigmaEst(y_expt);
    end
    if ~exist('debug','var') || isempty(debug)
        debug = false;
    end

    param0 = guess;
    
    kb_xind = find(mod(x_bp,1000)==0);
    x = x_bp(mean([kb_xind(1:end-1),kb_xind(2:end)],2));
    
    % Display off turns off the optimization complete statement. Lower and
    % upper bounds are empty matrices.
    opts = optimset('Display','off');
    [fitparam,~,~,~,~,~,Jacob] = lsqnonlin( @fitdiff, param0, [],[], opts);
    yfit = Theory(x_bp,fitparam,gene_CDHO);
    
    Info = full((Jacob'*Jacob)/(sigma^2));
    param_err = sqrt(diag(inv(Info)))';
    
    function dy = fitdiff( param )
        
        y_th = Theory(x_bp,param,gene_CDHO);
        dy = ( y_expt - y_th );
        
        if debug
            
            figure(1);
            clf;
        
            plot( x, y_expt, 'r.' );
            hold on;
            
            plot( x, y_th, 'b.-' );
            
            drawnow;
            pause(.01);
            
        end
    end
end

function y = Theory(x_bp,param,gene_CDHO)

    CD = param(1);
    HO = param(2);
    startheight = param(3);
    
    genes_L = gene_CDHO(x_bp<0);
    genes_R = gene_CDHO(x_bp>=0);
    
    slopes_L = genes_L;
    slopes_L(genes_L >= 0) = CD;
    slopes_L(genes_L == -1) = HO;
    
    slopes_R = genes_R;
    slopes_R(genes_R >= 0) = -CD;
    slopes_R(genes_R == -1) = -HO;
    
    slopes = [slopes_L;slopes_R];
    
    y_bp = cumsum(slopes)+startheight;
    
    kb_xind = find(mod(x_bp,1000)==0);
    yind = mean([kb_xind(1:end-1),kb_xind(2:end)],2);
    
    y = y_bp(yind);
    
end
