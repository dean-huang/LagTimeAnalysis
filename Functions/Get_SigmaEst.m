% Uses a nearest neighbors estimate to obtain standard deviation.

function sigma = Get_SigmaEst( data )
    
    sigma = sqrt(mean((data(2:end)-data(1:end-1)).^2,'omitnan'))/sqrt(2);

end