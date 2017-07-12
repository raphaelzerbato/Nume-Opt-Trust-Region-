function  [a] = ESN(x,mod,sigma,eps)
% if eps<-1 && eps>1
%     display -1<eps<1 is not respected
% else
    if (x-mod)< 0
        a = ((pi*2)^(-0.5)* exp(-(((x-mod)/sigma).^2)/(2*(1+eps)^2)))/sigma;
    else 
        a = ((pi*2)^(-0.5)* exp(-(((x-mod)/sigma).^2)/(2*(1-eps)^2)))/sigma;
    end
% end
end