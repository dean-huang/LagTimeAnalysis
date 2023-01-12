%% medcutfilt
% y: data to be filtered
% l window size for medfilt2
% fraction of extreme points to be removed.
function y = medcutfilt( y, l, p )

if sum(isnan(y)) == 0
    ny = numel( y );
else
    ny = sum(~isnan(y));
    y(isnan(y)) = inf;
end

y_mn = medfilt2( y, [l,l], 'symmetric' );

[~,ord] = sort( abs( y-y_mn ) );

inder = round((1-p)*ny);

if inder > ny
    inder = ny;
elseif inder<1
    inder = 1;
end


y(ord(inder:end)) = nan; 

end

