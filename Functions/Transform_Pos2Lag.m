function [tau] = Transform_Pos2Lag(x,logy,xpts,kG,debug)

kGmins = 60*kG;

[~,~,logyfit] = Fit_SingleSlope(x,logy,xpts,[],debug);

logyfitnorm = logyfit-max(logyfit);

tau = -kGmins^(-1)*logyfitnorm;

end

