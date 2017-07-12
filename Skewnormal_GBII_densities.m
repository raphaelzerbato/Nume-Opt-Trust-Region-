clear ''all'';
close all;
clc;

%% contexte
a = [0.50,0.50;0.50,0.50];

moyenne=[50,5];
ecart_type = [110,50];
qgauche = [0,0.1];
qdroite = [0.9];

mu=[moyenne(1,1),moyenne(1,2)];
sigma = [ecart_type(1,1),ecart_type(1,2)];
p= [qgauche(1,1),qgauche(1,2)];
q= [qdroite(1,1)];


muini=mu(1,2);
sigmaini =sigma(1,2);
pini=p(1,2);
qini= q;

pi1= [0.5,0.5];
o= xlsread('tablo time series.xlsx','TIME_SERIES', 'C2:C151');

one = [1;1];
N=2;
maxit = 150;
limit = 0.0001;
T = length(o);
gammma = ones(N,T);
digamma= ones(2*N,T-1);

%% definition du temps 
tps(1,1)=1;
for t = 1:T-1
    tps(t+1,1)=tps(t,1)+1;
end
% graf histo1
%histo
densite_essai=(@(x)(0.5* ESN(x,mu(1,1),sigma(1,1),p(1,1)))...
            +0.5*((mu(1,2)*x.^(mu(1,2)*p(1,2)-1))./(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(x/sigma(1,2)).^mu(1,2)).^(p(1,2)+q))));

densgb= densite_essai(o);
        
figure;
title('tentative de param')
hold on 
histogram((o),50,'Normalization','pdf');
scatter((o),(densgb))
axis([0 inf 0 0.0005])
hold off

%% 1 er histo
figure;
hold on
art = 25;
cc = zeros(T,3);
s_t = scatter(tps,o,art,cc,'filled');
colorbar
caxis([0 1])
drawnow
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pass 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha pass
f_1 = ESN(o(1,1),mu(1,1),sigma(1,1),p(1,1));
f_2 =((mu(1,2)*o(1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
b(1,1)= f_1;
b(2,2)= f_2; 


format long
a_1 = pi1*b;

%scaling alpha_1
format long
s = a_1*one;
salpha_1 = a_1./s;

%alpha pass1
format long
salphainter=zeros(N,T);
alpha(:,1)= a_1.';
salpha(:,1)= salpha_1.';

for k=2:T
    format long
    f_1 = ESN(o(k,1),mu(1,1),sigma(1,1),p(1,1));
    f_2 =((mu(1,2)*o(k,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(k,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
    
    b(1,1)= f_1;
    b(2,2)= f_2; 
    alpha(:,k) = salpha(:,k-1).'* a * b;
 
  %scaling alpha pass
    format long  
    s(1,k)= alpha(:,k).'*one;
    format long
    salphainter(:,k)= alpha(:,k)./s(1,k);
    
        if salphainter(1,k)<1e-135 
            salphainter(1,k)=0;
        elseif salphainter(2,k)<1e-135
            salphainter(2,k)=0;
        end
    salpha(:,k)= salphainter(:,k);
end

%% Beta pass1
format long
s_b = 1./ s;
format long
sbeta(:,T)=[s_b(1,T);s_b(1,T)];

for j = T-1:-1:1
    format long
    f_1 = ESN(o(j+1,1),mu(1,1),sigma(1,1),p(1,1));
    f_2 =((mu(1,2)*o(j+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(j+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
    b(1,1)= f_1;
    b(2,2)= f_2;
    format long
    beta(:,j)= a*b*sbeta(:,j+1);
  %beta scaling
    format long
    sbeta(:,j)= beta(:,j).*s_b(1,j);
end
format long
sbeta(1,T) = 1./s(1,T);

%% gamma1

for k = 1:(T-1)
    format long
    f_1 = ESN(o(k+1,1),mu(1,1),sigma(1,1),p(1,1));
    f_2 =((mu(1,2)*o(k+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(k+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
    
    b(1,1)= f_1;
    b(2,2)= f_2;  
    
    format long
    sa_1 = [salpha(1,k),0;0,salpha(1,k)];


    matb = vec2mat(sbeta(:,k+1),4);
    matb_1 = reshape(matb,2,2);
    matb_1(2,2)=matb_1(2,1);
    matb_1(2,1)= 0;

    tra = vec2mat(a(1,:),4);
    tra_1 = reshape (tra,2,2);
    tra_1(2,2) = tra_1(2,1);
    tra_1(2,1)=0;

    format long
    digam = sa_1*matb_1*tra_1*b*one;


    gammma(1,k)=(one.'*digam);
    
 
    sa_2 = [salpha(2,k),0;0,salpha(2,k)];

    tra2 = vec2mat(a(2,:),4);
    tra_2 = reshape (tra2,2,2);
    tra_2(2,2) = tra_2(2,1);
    tra_2(2,1)=0;

    format long
    digam2 =  sa_2*matb_1*tra_2*b*one;

    gammma(2,k)= (one.'*digam2);
       
    % reshaping digamma
    digam_1 = vec2mat(digam,4,[digam2(1,1),digam2(2,1)]);
    digamma(:,k) = reshape(digam_1,4,1);

end
% gamma_T
gammma(:,T)= alpha(:,T)./ s(1,T);
%% calcule P(O/lambda)
P(1,1) = -((log10(s))*ones(T,1))

%% restimation1

%pi
pi1= gammma(:,1).';

% A

for e = 1:2
    a(1,e)= (digamma(e,:)*ones(T-1,1)) / (gammma(1,(1:T-1))*ones(T-1,1));
    a(2,e)= (digamma(e+2,:)*ones(T-1,1))/(gammma(2,(1:T-1))*ones(T-1,1));
end


% B
%%%loi skn
ogamma1= [o,gammma(1,:)'];
S1=sortrows(ogamma1,1);

%%% k int 1 ere loi
G=zeros(100,T-1);
g=zeros(100,1);

for j = 1:T-1
    lp=linspace(S1(j,1),S1(j+1,1))';
    for k = 1:100
        g(k,1) = 1/4*(((ones(1,j)*((S1(1:j,1)-lp(k,1)).^2.*S1((1:j),2)))^(1/3))+((ones(1,T-j)*((S1(j+1:T,1)-lp(k,1)).^2.*S1((j+1:T),2)))^(1/3)))^3;
    end
    G(:,j)= g(:,1);
end
M =  min(G,[],1);
m = min(M);
[row1,column1] = find(G==m);
column1 =  column1(1,1);
f1minimod = S1(column1,1)+((S1(column1+1,1)-S1(column1,1))*(column1/100));

numepsilone =((ones(1,column1)*((((S1(1:column1,1)-f1minimod)).^2).*(S1(1:column1,2))))^(1/3))-((ones(1,T-column1)*((((S1(column1+1:T,1)-f1minimod)).^2).*(S1(column1+1:T,2))))^(1/3));
denepsilone = ((ones(1,column1)*((((S1(1:column1,1)-f1minimod)).^2).*(S1(1:column1,2))))^(1/3))+((ones(1,T-column1)*((((S1(column1+1:T,1)-f1minimod)).^2).*(S1(column1+1:T,2))))^(1/3));
f1epsesti =  numepsilone / denepsilone;

f1sigmaesti = ((1/(4*(gammma(1,(1:T))*ones(T,1))))*(denepsilone)^3);
MLE1(1,2)= -((gammma(1,(1:T))*ones(T,1))/2)*(1+log10(f1sigmaesti^2));

mu(1,1)= f1minimod;
sigma(1,1)=(f1sigmaesti)^(1/2);
p(1,1)= f1epsesti;


%%% loi GB2
f2=@(x)-sum(sum(gammma(2,:),2)*(log(x(1))-log(x(2))-log(gamma(x(3)))-log(gamma(x(4)))+log(gamma(x(3)+x(4))))+gammma(2,:)*(x(1)*x(3)-1)*log(o/x(2))-gammma(2,:)*(x(3)+x(4))*log(1+(o/x(2)).^x(1)));
%options = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','bfgs','MaxFunEvals',700);
lb = zeros(1,4);
ub = +Inf;
xfinal2 = fmincon(f2,[muini,sigmaini,pini,qini],[],[],[],[],lb,ub)




mu(1,2) = xfinal2(1,1);
sigma(1,2) = (xfinal2(1,2));
p(1,2) = xfinal2(1,3);
q = xfinal2(1,4);

%%graph%%
figure_2 = figure;
subplot(2,2,3);
hold on
art = 25;
cc = zeros(T,3);
gammabis2 = gammma';
onne = gammabis2(:,1);
twoo = gammabis2(:,2);
res = onne - twoo;

for temp=1:T
    if res(temp,1) > 0
        cc(temp,:)= [0 0 (1+res(temp))/2 ];
    else
        cc(temp,:)= [0 (1+res(temp))/2  0];
    end
end

s_t = scatter(tps,log(o),art,cc,'filled');
%colorbar
%caxis([0 1])
ylabel('Log(data)')
xlabel('Temps')
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
set(gca,'XTick',[]);
caxis([0 1])
drawnow
title('Observations: première itération','FontSize',8)
caxis([0 1])
drawnow
hold off
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 ème pass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha pass 2
f_1 = ESN(o(1,1),mu(1,1),sigma(1,1),p(1,1));
f_2 =((mu(1,2)*o(1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
b(1,1)= f_1;
b(2,2)= f_2; 


format long
a_1 = pi1*b;

%scaling alpha_2
format long
s = a_1*one;
salpha_1 = a_1./s;

%alpha pass 2
format long
alpha(:,1)= a_1.';
salpha(:,1)= salpha_1.';

for k=2:T
    format long
    f_1 = ESN(o(k,1),mu(1,1),sigma(1,1),p(1,1));
    f_2 =((mu(1,2)*o(k,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(k,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
    
    b(1,1)= f_1;
    b(2,2)= f_2; 
    alpha(:,k) = salpha(:,k-1).'* a * b;
 
  %scaling alpha pass
    format long  
    s(1,k)= alpha(:,k).'*one;
    format long
    salphainter(:,k)= alpha(:,k)./s(1,k);
    
         if salphainter(1,k)<1e-135
            salphainter(1,k)=0;
        elseif salphainter(2,k)<1e-135
            salphainter(2,k)=0;
        end
    salpha(:,k)= salphainter(:,k);
end

%% Beta pass 2 
format long
s_b = 1./ s;
format long
sbeta(:,T)=[s_b(1,T);s_b(1,T)];

for j = T-1:-1:1
    format long
    f_1 = ESN(o(j+1,1),mu(1,1),sigma(1,1),p(1,1));
    f_2 =((mu(1,2)*o(j+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(j+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
    b(1,1)= f_1;
    b(2,2)= f_2;
    format long
    beta(:,j)= a*b*sbeta(:,j+1);
  %beta scaling
    format long
    sbeta(:,j)= beta(:,j).*s_b(1,j);
end
format long
sbeta(1,T) = 1./s(1,T);

%% gamma 2 

for k = 1:(T-1)
    format long
    f_1 = ESN(o(k+1,1),mu(1,1),sigma(1,1),p(1,1));
    f_2 =((mu(1,2)*o(k+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(k+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
    
    b(1,1)= f_1;
    b(2,2)= f_2;  
    
    format long
    sa_1 = [salpha(1,k),0;0,salpha(1,k)];


    matb = vec2mat(sbeta(:,k+1),4);
    matb_1 = reshape(matb,2,2);
    matb_1(2,2)=matb_1(2,1);
    matb_1(2,1)= 0;

    tra = vec2mat(a(1,:),4);
    tra_1 = reshape (tra,2,2);
    tra_1(2,2) = tra_1(2,1);
    tra_1(2,1)=0;

    format long
    digam = sa_1*matb_1*tra_1*b*one;


    gammma(1,k)=(one.'*digam);
    
 
    sa_2 = [salpha(2,k),0;0,salpha(2,k)];

    tra2 = vec2mat(a(2,:),4);
    tra_2 = reshape (tra2,2,2);
    tra_2(2,2) = tra_2(2,1);
    tra_2(2,1)=0;

    format long
    digam2 =  sa_2*matb_1*tra_2*b*one;

    gammma(2,k)= (one.'*digam2);
       
    % reshaping digamma
    digam_1 = vec2mat(digam,4,[digam2(1,1),digam2(2,1)]);
    digamma(:,k) = reshape(digam_1,4,1);

end
% gamma_T
gammma(:,T)= alpha(:,T)./ s(1,T);
%% calcule P(O/lambda)
P(1,2) = -((log10(s))*ones(T,1))

%% restimation 2

%pi
pi1= gammma(:,1).';

% A

for e = 1:2
    a(1,e)= (digamma(e,:)*ones(T-1,1)) / (gammma(1,(1:T-1))*ones(T-1,1));
    a(2,e)= (digamma(e+2,:)*ones(T-1,1))/(gammma(2,(1:T-1))*ones(T-1,1));
end


% B
%%%loi skn
ogamma1= [o,gammma(1,:)'];
S1=sortrows(ogamma1,1);

%%% k int 1 ere loi
G=zeros(100,T-1);
g=zeros(100,1);

for j = 1:T-1
    lp=linspace(S1(j,1),S1(j+1,1))';
    for k = 1:100
        g(k,1) = 1/4*(((ones(1,j)*((S1(1:j,1)-lp(k,1)).^2.*S1((1:j),2)))^(1/3))+((ones(1,T-j)*((S1(j+1:T,1)-lp(k,1)).^2.*S1((j+1:T),2)))^(1/3)))^3;
    end
    G(:,j)= g(:,1);
end
M =  min(G,[],1);
m = min(M);
[row1,column1] = find(G==m);
column1 =  column1(1,1);
f1minimod = S1(column1,1)+((S1(column1+1,1)-S1(column1,1))*(column1/100));

numepsilone =((ones(1,column1)*((((S1(1:column1,1)-f1minimod)).^2).*(S1(1:column1,2))))^(1/3))-((ones(1,T-column1)*((((S1(column1+1:T,1)-f1minimod)).^2).*(S1(column1+1:T,2))))^(1/3));
denepsilone = ((ones(1,column1)*((((S1(1:column1,1)-f1minimod)).^2).*(S1(1:column1,2))))^(1/3))+((ones(1,T-column1)*((((S1(column1+1:T,1)-f1minimod)).^2).*(S1(column1+1:T,2))))^(1/3));
f1epsesti =  numepsilone / denepsilone;

f1sigmaesti = ((1/(4*(gammma(1,(1:T))*ones(T,1))))*(denepsilone)^3);
MLE1(1,2)= -((gammma(1,(1:T))*ones(T,1))/2)*(1+log10(f1sigmaesti^2));

mu(1,1)= f1minimod;
sigma(1,1)=(f1sigmaesti)^(1/2);
p(1,1)= f1epsesti;


%%% loi GB2
f2=@(x)-sum(sum(gammma(2,:),2)*(log(x(1))-log(x(2))-log(gamma(x(3)))-log(gamma(x(4)))+log(gamma(x(3)+x(4))))+gammma(2,:)*(x(1)*x(3)-1)*log(o/x(2))-gammma(2,:)*(x(3)+x(4))*log(1+(o/x(2)).^x(1)));
%options = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','bfgs','MaxFunEvals',700);
lb = zeros(1,4);
ub = +Inf;
xfinal2 = fmincon(f2,[muini,sigmaini,pini,qini],[],[],[],[],lb,ub)

mu(1,2) = xfinal2(1,1);
sigma(1,2) = (xfinal2(1,2));
p(1,2) = xfinal2(1,3);
q = xfinal2(1,4);

%%graph%%
figure;
         hold on
         art = 25;
         cc = zeros(T,3);
         gammabis2 = gammma';
         onne = gammabis2(:,1);
         twoo = gammabis2(:,2);
         res = onne - twoo;
        
         for temp=1:T
             if res(temp,1) > 0
                 cc(temp,:)= [0 0 (1+res(temp))/2 ];
             else
                 cc(temp,:)= [0 (1+res(temp))/2  0];
             end
         end

         s_t = scatter(tps,log(o),art,cc,'filled');
         colorbar
         caxis([0 1])
         drawnow
         hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 ème pass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leme pass
for L= 2:maxit
    %% alpha pass
    f_1 = ESN(o(1,1),mu(1,1),sigma(1,1),p(1,1));
    f_2 =((mu(1,2)*o(1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
    b(1,1)= f_1;
    b(2,2)= f_2;

    
    format long
    a_1 = pi1*b;
    
    %scaling alpha_2
    format long
    s = a_1*one;
    salpha_1 = a_1./s;
    
    %alpha pass 2
    format long
    alpha(:,1)= a_1.';
    salpha(:,1)= salpha_1.';
    
    for k=2:T
        format long
        f_1 = ESN(o(k,1),mu(1,1),sigma(1,1),p(1,1));
        f_2 =((mu(1,2)*o(k,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(k,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
        
        b(1,1)= f_1;
        b(2,2)= f_2;
        alpha(:,k) = salpha(:,k-1).'* a * b;
        
        %scaling alpha pass
        format long
        s(1,k)= alpha(:,k).'*one;
        format long
        salphainter(:,k)= alpha(:,k)./s(1,k);
    
        if salphainter(1,k)<1e-135 
            salphainter(1,k)=0;
        elseif salphainter(2,k)<1e-135
            salphainter(2,k)=0;
        end
    salpha(:,k)= salphainter(:,k);
    
    end
    
    %% Beta pass 2
    format long
    s_b = 1./ s;
    format long
    sbeta(:,T)=[s_b(1,T);s_b(1,T)];
    
    for j = T-1:-1:1
        format long
        f_1 = ESN(o(j+1,1),mu(1,1),sigma(1,1),p(1,1));
        f_2 =((mu(1,2)*o(j+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(j+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
        b(1,1)= f_1;
        b(2,2)= f_2;
        format long
        beta(:,j)= a*b*sbeta(:,j+1);
        %beta scaling
        format long
        sbeta(:,j)= beta(:,j).*s_b(1,j);
    end
    format long
    sbeta(1,T) = 1./s(1,T);
    
    %% gamma 2
    
    for k = 1:(T-1)
        format long
        f_1 = ESN(o(k+1,1),mu(1,1),sigma(1,1),p(1,1));
        f_2 =((mu(1,2)*o(k+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q)/gamma(p(1,2)+q))*(1+(o(k+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q)));
        
        b(1,1)= f_1;
        b(2,2)= f_2;
        
        format long
        sa_1 = [salpha(1,k),0;0,salpha(1,k)];
        
        
        matb = vec2mat(sbeta(:,k+1),4);
        matb_1 = reshape(matb,2,2);
        matb_1(2,2)=matb_1(2,1);
        matb_1(2,1)= 0;
        
        tra = vec2mat(a(1,:),4);
        tra_1 = reshape (tra,2,2);
        tra_1(2,2) = tra_1(2,1);
        tra_1(2,1)=0;
        
        format long
        digam = sa_1*matb_1*tra_1*b*one;
        
        
        gammma(1,k)=(one.'*digam);
        
        
        sa_2 = [salpha(2,k),0;0,salpha(2,k)];
        
        tra2 = vec2mat(a(2,:),4);
        tra_2 = reshape (tra2,2,2);
        tra_2(2,2) = tra_2(2,1);
        tra_2(2,1)=0;
        
        format long
        digam2 =  sa_2*matb_1*tra_2*b*one;
        
        gammma(2,k)= (one.'*digam2);
        
        % reshaping digamma
        digam_1 = vec2mat(digam,4,[digam2(1,1),digam2(2,1)]);
        digamma(:,k) = reshape(digam_1,4,1);
        
    end
    % gamma_T
    gammma(:,T)= alpha(:,T)./ s(1,T);
    %% calcule P(O/lambda)
    P(1,L+1) = -((log10(s))*ones(T,1))    
    %% restimation 2
    
    %pi
    pi1= gammma(:,1).';
    
    % A
    
    for e = 1:2
        a(1,e)= (digamma(e,:)*ones(T-1,1)) / (gammma(1,(1:T-1))*ones(T-1,1));
        a(2,e)= (digamma(e+2,:)*ones(T-1,1))/(gammma(2,(1:T-1))*ones(T-1,1));
    end
    
    
    % B
    %%%loi skn
    ogamma1= [o,gammma(1,:)'];
    S1=sortrows(ogamma1,1);
    
    %%% k int 1 ere loi
    G=zeros(100,T-1);
    g=zeros(100,1);
    
    for j = 1:T-1
        lp=linspace(S1(j,1),S1(j+1,1))';
        for k = 1:100
            g(k,1) = 1/4*(((ones(1,j)*((S1(1:j,1)-lp(k,1)).^2.*S1((1:j),2)))^(1/3))+((ones(1,T-j)*((S1(j+1:T,1)-lp(k,1)).^2.*S1((j+1:T),2)))^(1/3)))^3;
        end
        G(:,j)= g(:,1);
    end
    M =  min(G,[],1);
    m = min(M);
    [row1,column1] = find(G==m);
    column1 =  column1(1,1);
    f1minimod = S1(column1,1)+((S1(column1+1,1)-S1(column1,1))*(column1/100));
    
    numepsilone =((ones(1,column1)*((((S1(1:column1,1)-f1minimod)).^2).*(S1(1:column1,2))))^(1/3))-((ones(1,T-column1)*((((S1(column1+1:T,1)-f1minimod)).^2).*(S1(column1+1:T,2))))^(1/3));
    denepsilone = ((ones(1,column1)*((((S1(1:column1,1)-f1minimod)).^2).*(S1(1:column1,2))))^(1/3))+((ones(1,T-column1)*((((S1(column1+1:T,1)-f1minimod)).^2).*(S1(column1+1:T,2))))^(1/3));
    f1epsesti =  numepsilone / denepsilone;
    
    f1sigmaesti = ((1/(4*(gammma(1,(1:T))*ones(T,1))))*(denepsilone)^3);
    MLE1(1,2)= -((gammma(1,(1:T))*ones(T,1))/2)*(1+log10(f1sigmaesti^2));
    
    mu(1,1)= f1minimod;
    sigma(1,1)=(f1sigmaesti)^(1/2);
    p(1,1)= f1epsesti;
    
    
    %%% loi GB2
    f2=@(x)-sum(sum(gammma(2,:),2)*(log(x(1))-log(x(2))-log(gamma(x(3)))-log(gamma(x(4)))+log(gamma(x(3)+x(4))))+gammma(2,:)*(x(1)*x(3)-1)*log(o/x(2))-gammma(2,:)*(x(3)+x(4))*log(1+(o/x(2)).^x(1)));
    lb = zeros(1,4);
    ub = +Inf;
    xfinal2 = fmincon(f2,[muini,sigmaini,pini,qini],[],[],[],[],lb,ub)
    
    mu(1,2) = xfinal2(1,1);
    sigma(1,2) = (xfinal2(1,2));
    p(1,2) = xfinal2(1,3);
    q = xfinal2(1,4);
    
     
     %% graphique
     if L~=maxit
         if mod(L,10)==0
             figure;
             hold on
             art = 25;
             cc = zeros(T,3);
             gammabis2 = gammma';
             onne = gammabis2(:,1);
             twoo = gammabis2(:,2);
             res = onne - twoo;
             for temp=1:T
                 if res(temp,1) > 0
                     cc(temp,:)= [0 0 (1+res(temp))/2 ];
                 else
                     cc(temp,:)= [0 (1+res(temp))/2  0];
                 end
             end
             
             s_t = scatter(tps,log(o),art,cc,'filled');
             colorbar
             caxis([0 1])
             drawnow
             hold off
         end
     elseif L==maxit 
        figure(figure_2);
            subplot(2,2,4) 
            hold on
             art = 25;
             cc = zeros(T,3);
             gammabis2 = gammma';
             onne = gammabis2(:,1);
             twoo = gammabis2(:,2);
             res = onne - twoo;
             for temp=1:T
                 if res(temp,1) > 0
                     cc(temp,:)= [0 0 (1+res(temp))/2 ];
                 else
                     cc(temp,:)= [0 (1+res(temp))/2  0];
                 end
             end
             s_t = scatter(tps,log(o),art,cc,'filled');
             %colorbar
             set(gca, 'FontName', 'Arial')
             set(gca, 'FontSize', 8)
             ylabel('Log(data)')
             xlabel('Temps')
             set(gca,'XTick',[]);
             caxis([0 1])
             drawnow
             title('Observations: dernière itération','FontSize',8)
             hold off
     end
end
%% graph de fin

% Taille
taille = [(gammma(1,:)*ones(T,1)/length(o)); ((gammma(2,:)*ones(T,1))/length(o))];

%graf de la convergences des probas
figure_4 = figure
hold on
iteration =  ones(1,length(P));
for l = 1 : length(P)-1
    iteration(1,l+1)=iteration(1,l)+1;
end    
proba = plot(iteration,P(1,:))
drawnow
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('P(O/lambda)')
xlabel('Iterations')
title('convergence des proba','FontSize',8)
hold off

%% graf de verification 

%densité gb2
a2= xfinal2(1,1);
b2 =xfinal2(1,2);
p2= xfinal2(1,3);
q2= xfinal2(1,4);

funct1 =@(z)(taille(1,1)*ESN(z,mu(1,1),sigma(1,1),p(1,1))...
            +taille(2,1)*((a2*z.^(a2*p2-1))./(b2^(a2*p2)*(gamma(p2)*gamma(q2)/gamma(p2+q2))*(1+(z/b2).^a2).^(p2+q2))));
densgb = funct1(o);

figure_1 = figure;
hold on 
subplot(2,2,1)
histogram((o),150,'Normalization','pdf','FaceColor','b')
hold on
h=scatter((o),densgb,5,'filled','m');
axis([-inf,inf,0,0.01]);
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('density')
xlabel('observation en million')
hold off
title('Donnée et densité', 'FontSize', 8)
hold off



%graf des densites
figure(figure_2)
subplot(2,2,1)
hold on
ezplot(@(z)(taille(1,1)*ESN(z,mu(1,1),sigma(1,1),p(1,1))...
            +taille(2,1)*((a2*z.^(a2*p2-1))./(b2^(a2*p2)*(gamma(p2)*gamma(q2)/gamma(p2+q2))*(1+(z/b2).^a2).^(p2+q2)))),[0,1000,0,0.01]);
title('Fonction de densité', 'FontSize', 8)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('density')
xlabel('observation en million')
hold off


%% influence des PDF Les deux CDF sont comparés   
%%% construciton de la CDF GB2 

y=(o./sigma(1,2)).^mu(1,2);
z = y./(1+y);
GB2cdf = betainc(z,p(1,2),q);


%%% construction de la CDF Skewn
vect1=zeros(T,1);
for i=1:T
    vect1(i,1) = SDESNCDF(o(i,1),mu(1,1),sigma(1,1),p(1,1));
end
skncdf = vect1;


testpdf=zeros(T,1);
for i=1:T
    testpdf(i,1) = taille(2,1)*GB2cdf(i,1)+taille(1,1)*skncdf(i,1);
end

figure(figure_2)
subplot(2,2,2)
hold on
scatter(o,testpdf,10,'o','filled')
scatter(o,vect1,'.')
scatter(o,GB2cdf,5,'d')
drawnow
axis([0,max(o),0,1]);
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('probability')
xlabel('observation en million')
legend('CdF tot','CdF 1','CdF 2')
title('CDF', 'FontSize', 8)
hold off
print('figure_2','-djpeg')
print ('figure_2','-dpdf','-r600')
%% comparaison PDF et test  

parampdf =@(vari) (taille(1,1)*(SDESNCDF(vari,mu(1,1),sigma(1,1),p(1,1)))...
    + taille(2,1)*GB2CDF(vari,mu(1,2),sigma(1,2),p(1,2),q));

[empdist,x] = ecdf(o);


figure(figure_1)
subplot(2,2,2)
hold on
ezplot(parampdf,[0,max(o),0,1])
plot(x,empdist)
title('CDF empirique et CDF modèle', 'FontSize',8)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('probability')
xlabel('observation en million')
legend('CdF param','CdF empriq')
hold off

%% test des pseudo residual 

unipseudo=zeros(T,1);
for i = 1:T
    unipseudo(i,1)= parampdf(o(i,1));   
end
npseudores = zeros(T,1);
for i=1:T
    npseudores(i,1) = norminv(unipseudo(i,1),0,1);
end
%%% normal QQ plot
Quantiltest1 = [0.10:0.01:0.95]';
zq=zeros(length(Quantiltest1),1);

for i=1:length(Quantiltest1)
    zq(i,1)=norminv(Quantiltest1(i,1),0,1);
end
zquant = quantile(npseudores,Quantiltest1);

figure(figure_1);
subplot(2,2,4)
qqplot(zq,zquant)
title('QQ plot pseudoresi', 'FontSize',8)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
xlabel('Std Normal Quantile')
ylabel('z quantile')
hold off

[hchi,pchi,stat] = chi2gof(zquant,'Alpha',0.05); % test chi

[hswt,pswtest]=swtest(zquant,0.05); % test swtest

[httest,pttest] = ttest(zquant,0.05); % Ttest 

AIC = -2*(-P(1,L))+20;
BIC = -2*(-P(1,L))+10*log(T);



format short % tableau des resultats des test
parametre1 = {'0 -> Normal accepted at 5%';'P values'};
Chi_square_test = [hchi;pchi];
Shapiro_Wilk_test = [hswt;pswtest];
table_parametres1 = table(Chi_square_test,Shapiro_Wilk_test,'RowNames',parametre1)

figure(figure_4);
t = uitable(figure_4,'Data',{'Chi_square_test',hchi,pchi;...
    'Shapiro_Wilk_test',hswt,pswtest;'Ttest',httest,pttest},'Position',[270 320 258 70]);
t.ColumnName = {'test','H1 true 5%','P values'};
t.ColumnEditable = true;
abt = uitable(figure_4,'Data',{AIC,BIC},'Position',[270 230 200 70]);
abt.ColumnName = {'AIC criteria','BIC criteria'};
abt.ColumnEditable= true;
hold off
print('figure_4','-djpeg')
print ('figure_4','-dpdf','-r600')




%% creation des QQ plots pour les donnees 
Quantiltest = [0.01:0.05:0.99]';
EmQ=zeros(length(Quantiltest),1);

for i=1:length(Quantiltest)
    EmQ(i,1)= (1/taille(1,1))*(generatorESN(Quantiltest(i,1),mu(1,1),sigma(1,1),p(1,1)))...
      +(1/taille(2,1))*(GB2qtl(Quantiltest(i,1),mu(1,2),sigma(1,2),p(1,2),q));
end

quantdata = quantile(o,Quantiltest);

figure(figure_1)
subplot(2,2,3)
qqplot(EmQ,quantdata);
title('QQ plot', 'FontSize',8)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
xlabel('Theoretical Quantile')
ylabel('Data quantile')
hold off
print('figure_1','-djpeg')
print ('figure_1','-dpdf','-r600')
%  %%
% 
% Quantiltest = [0.01:0.05:0.99]';
% 
% 
% EmQ1=zeros(length(Quantiltest),1);
% for i=1:length(Quantiltest)
%     EmQ1(i,1)= (1/taille(1,1))*(1/2+norminv(Quantiltest(i,1),mu(1,1),sigma(1,1)))...
%         +(1/taille(2,1))*(1/2+GB2qtl(Quantiltest(i,1),mu(1,2),sigma(1,2),p(1,2),q));
% end
% 
% quantdata1 = quantile(o,Quantiltest);
% 
% figure
% qqplot(quantdata1,EmQ1);
% title('QQ plot', 'FontSize',8)
% set(gca, 'FontName', 'Arial')
% set(gca, 'FontSize', 8)
% ylabel('Theoretical Quantile')
% xlabel('Data quantile')
% hold off

%% creation de la figure des gammabis
figure_3 = figure
%nt = datetime(1999,01,01) + calmonths(0:149);
t =(1:1:150);
%hold on
%plot(t,gammabis2(:,1),'--k')
hold on 
plot(t,gammabis2(:,2),'--k')
hold on
ar1 = area([27 35], [1 1]); % creation of crisis shades 1
set(ar1, 'FaceColor', [.5 .5 .5]);
set(ar1,'EdgeColor', 'none')
drawnow; 
pause(0.05);  % This needs to be done for transparency to work
ar1.Face.ColorType = 'truecoloralpha';
ar1.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3

br1 = area([108 126], [1 1]); % creation of crisis shades 1
set(br1, 'FaceColor', [.5 .5 .5]);
set(br1,'EdgeColor', 'none')
drawnow; 
pause(0.05);  % This needs to be done for transparency to work
br1.Face.ColorType = 'truecoloralpha';
br1.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3
legend('loi 2')
title('régime', 'FontSize',10)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('Probabilité état')
xlabel('Temps')
hold off
print('figure_3','-djpeg')
print ('figure_3','-dpdf','-r600')

%% creation table de resume des parametres
format short
parametre = {'Loi1final';'Loi2final';'Loi1ini';'Loi2ini'};
A = [mu(1,1);mu(1,2);moyenne(1,1);moyenne(1,2)];
B = [sigma(1,1);sigma(1,2);ecart_type(1,1);ecart_type(1,2)];
P = [p(1,1);p(1,2);qgauche(1,1);qgauche(1,2)];
Q = [0;q;0;qdroite];
table_parametres = table(A,B,P,Q,'RowNames',parametre)

figure_5 = figure;
format short
t = uitable(figure_5,'Data',{'Loi1final',mu(1,1),sigma(1,1),p(1,1),'NA'...
    ;'loi2final',mu(1,2),sigma(1,2),p(1,2),q;'loi1ini',moyenne(1,1),...
    ecart_type(1,1),qgauche(1,1),'NA';'loi2ini',moyenne(1,2),...
    ecart_type(1,2),qgauche(1,2),qdroite},'Position',[85 250 410 102]);
t.ColumnName = {'loi','A','B','P/eps','Q'};
t.ColumnEditable = true;
ax = gca;
set(ax, 'Visible','off')
hold off
%% creation des tablo crise non crise
pbnocriseL1=(sum(gammabis2(1:26,1))+sum(gammabis2(35:107,1))...
    +sum(gammabis2(127:150,1)))/(150-27);

pbnocriseL2=(sum(gammabis2(1:26,2))+sum(gammabis2(35:107,2))...
    +sum(gammabis2(127:150,2)))/(150-27);

pbcriseL1=(sum(gammabis2(27:34,1))+sum(gammabis2(108:126,1)))/27;

pbcriseL2=(sum(gammabis2(27:34,2))+sum(gammabis2(108:126,2)))/27;

format short
tete = {'Crise';'Non-Crise';'Taille totale'};
Loi_1 = [pbcriseL1;pbnocriseL1;taille(1,1)];
Loi_2 = [pbcriseL2;pbnocriseL2;taille(2,1)];
table_crise = table(Loi_1,Loi_2,'RowNames',tete)


figure(figure_5);
format short
cris = uitable(figure_5,'Data',{'Loi_1',pbcriseL1,pbnocriseL1,taille(1,1)...
    ;'Loi_2',pbcriseL2,pbnocriseL2,taille(2,1)},'Position',[85 120 335 75]);
cris.ColumnName = {'Loi','Crise','Non Crise','Taille totale'};
cris.ColumnEditable = true;
ax = gca;
set(ax, 'Visible','off')
hold off
print('figure_5','-djpeg')
print ('figure_5','-dpdf','-r600')
mydoc;
