function a = GB2CDF(o,mod,std,p,q)
y=(o./std).^ mod;
x = y/(1+y);
a = betainc(x,p,q);
end