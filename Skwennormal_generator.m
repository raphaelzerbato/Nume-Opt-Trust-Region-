function res = generatorESN(p,mod,sigma,eps)
if p >= 0 && p < ((1+eps)/2)
    po = mod+sigma*((1+eps)*norminv((p/(1+eps)),0,1));
elseif p > ((1+eps)/2) &&  p < 1
    po = mod+sigma*((1-eps)*norminv(((p-eps)/(1-eps)),0,1));
end
res = po;
end