function [a] = SDESNCDF(y,mu,sttd,eps)
    x = (y-mu)/sttd;
    if x < 0
        a = (1+eps)*normcdf((x/(1+eps)),0,1);
    else 
        a = eps+(1-eps)*normcdf((x/(1-eps)),0,1);
    end

end
