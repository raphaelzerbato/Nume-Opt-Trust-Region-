function a = GB2qtl(o,mod,std,p,q)
lol=betainv(o,p,q);
y=lol/(1-lol);
a=std*y^(1/mod);
end