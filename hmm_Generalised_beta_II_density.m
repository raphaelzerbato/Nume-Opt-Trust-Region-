clear ''all'';
close all;
clc;

%% Mise en place
%%%%% Paramètres initiaux
% Le modèle estime deux régimes distinct utilisant des lois d'émission
% Beta généralisée de Type II à chaque régime.
% 'ar1' est la matrice de transition initiale donnée 
% 'A,B,P,Q' sont les vecteurs de paramètres initiaux choisit
%'pi'est la probabilité initiale régime 1 ou régime 2 
% maxit est le nombre d'itérations souhaité pour le modèle
% o est le vecteur de données desquelles sont estimés les régimes

ar1 = [0.50,0.50;0.50,0.50];

A=[2,10];
B = [10,1];
P = [0.05,0.09];
Q = [0.05,0.01];


mu=[A(1,1),A(1,2)];
sigma = [B(1,1),B(1,2)];
p= [P(1,1),P(1,2)];
q= [Q(1,1),Q(1,2)];
muini= mu;
sigmaini = sigma;
pini= p;
qini=q;

pi1= [0.5,0.5];
o= xlsread('tablo time series.xlsx','TIME_SERIES','C2:C151');
% date = xlsread(C:\Users\inkyiagan\Desktop\code mathlab1\MODEL fonctionel\GB2'tablo 1.xlsm','Fréquence des deux correction','A2:A151');
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
densite_essai=(@(x)(0.5*((mu(1,1)*x.^(mu(1,1)*p(1,1)-1))./(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(x/sigma(1,1)).^mu(1,1)).^(p(1,1)+q(1,1))))...
            +0.5*((mu(1,2)*x.^(mu(1,2)*p(1,2)-1))./(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(x/sigma(1,2)).^mu(1,2)).^(p(1,2)+q(1,2))))));

densgb= densite_essai(o);
        
figure_6=figure;
hold on 
histogram((o),150,'Normalization','pdf');
%scatter((o),(densgb))
title('histograme somme pertes mensuelles')
axis([0 inf 0 0.006])
ylabel('Densité')
xlabel('Observation')
hold off
print('figure_6','-djpeg')
print ('figure_6','-dpdf','-r600')
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

%% 1 er pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pass 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha pass
f_1 =((mu(1,1)*o(1,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(o(1,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
f_2 =((mu(1,2)*o(1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
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
alpha(:,1)= a_1.';
salpha(:,1)= salpha_1.';
salphainter = zeros(N,T);
for k=2:T
    format long
    f_1 =((mu(1,1)*o(k,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/...
        gamma(p(1,1)+q(1,1)))*(1+(o(k,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
    f_2 =((mu(1,2)*o(k,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(k,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
    b(1,1)= f_1;
    b(2,2)= f_2; 
    alpha(:,k) = salpha(:,k-1).'* ar1 * b;
 
  %scaling alpha pass
    format long  
    s(1,k)= alpha(:,k).'*one;
    format long
    salphainter(:,k)= alpha(:,k)./s(1,k);
    
    if salphainter(1,k)<1e-135 || isnan(salphainter(1,k))
        salphainter(1,k)=0;
    elseif salphainter(2,k)<1e-135|| isnan(salphainter(2,k))
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
    f_1 =((mu(1,1)*o(j+1,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(o(j+1,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
    f_2 =((mu(1,2)*o(j+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(j+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
    b(1,1)= f_1;
    b(2,2)= f_2;
    format long
    beta(:,j)= ar1*b*sbeta(:,j+1);
  %beta scaling
    format long
    sbeta(:,j)= beta(:,j).*s_b(1,j);
end
format long
sbeta(1,T) = 1./s(1,T);

%% gamma1

for k = 1:(T-1)
    format long
    f_1 =((mu(1,1)*o(k+1,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(o(k+1,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
    f_2 =((mu(1,2)*o(k+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(k+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
    
    b(1,1)= f_1;
    b(2,2)= f_2;  
    
    format long
    sa_1 = [salpha(1,k),0;0,salpha(1,k)];


    matb = vec2mat(sbeta(:,k+1),4);
    matb_1 = reshape(matb,2,2);
    matb_1(2,2)=matb_1(2,1);
    matb_1(2,1)= 0;

    tra = vec2mat(ar1(1,:),4);
    tra_1 = reshape (tra,2,2);
    tra_1(2,2) = tra_1(2,1);
    tra_1(2,1)=0;

    format long
    digam = sa_1*matb_1*tra_1*b*one;


    gammma(1,k)=(one.'*digam);
    
 
    sa_2 = [salpha(2,k),0;0,salpha(2,k)];

    tra2 = vec2mat(ar1(2,:),4);
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
Pb(1,1) = -((log10(s))*ones(T,1))
%% restimation1

%pi
pi1= gammma(:,1).';

% A

for e = 1:2
    ar1(1,e)= (digamma(e,:)*ones(T-1,1)) / (gammma(1,(1:T-1))*ones(T-1,1));
    ar1(2,e)= (digamma(e+2,:)*ones(T-1,1))/(gammma(2,(1:T-1))*ones(T-1,1));
end


% B
f1=@(x)-sum(sum(gammma(1,:),2)*(log(x(1))-log(x(2))-log(gamma(x(3)))-log(gamma(x(4)))+log(gamma(x(3)+x(4))))+gammma(1,:)*(x(1)*x(3)-1)*log(o/x(2))-gammma(1,:)*(x(3)+x(4))*log(1+(o/x(2)).^x(1)));
lb = zeros(1,4);
ub = +Inf;
xfinal1 = fmincon(f1,[muini(1,1),sigmaini(1,1),pini(1,1),qini(1,1)],[],[],[],[],lb,ub)



f2=@(x)-sum(sum(gammma(2,:),2)*(log(x(1))-log(x(2))-log(gamma(x(3)))-log(gamma(x(4)))+log(gamma(x(3)+x(4))))+gammma(2,:)*(x(1)*x(3)-1)*log(o/x(2))-gammma(2,:)*(x(3)+x(4))*log(1+(o/x(2)).^x(1)));
lb = zeros(1,4);
ub = +Inf;
xfinal2 = fmincon(f2,[muini(1,2),sigmaini(1,2),pini(1,2),qini(1,2)],[],[],[],[],lb,ub)

mu = [xfinal1(1,1),xfinal2(1,1)];
sigma = [(xfinal1(1,2)),(xfinal2(1,2))];
p = [xfinal1(1,3),xfinal2(1,3)];
q = [xfinal1(1,4),xfinal2(1,4)];

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
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Grande boucle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% L eme pass
for L= 2:maxit
    % alpha pass
    f_1 =((mu(1,1)*o(1,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(o(1,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
    f_2 =((mu(1,2)*o(1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
    b(1,1)= f_1;
    b(2,2)= f_2;
    
    
    format long
    a_1 = pi1*b;
    
    %scaling alpha
    format long
    s = a_1*one;
    salpha_1 = a_1./s;
    
    %alpha pass
    format long
    alpha(:,1)= a_1.';
    salpha(:,1)= salpha_1.';
    salphainter = zeros(N,T);
    for k=2:T
        format long
        f_1 =((mu(1,1)*o(k,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(o(k,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
        f_2 =((mu(1,2)*o(k,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(k,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
        b(1,1)= f_1;
        b(2,2)= f_2;
        alpha(:,k) = salpha(:,k-1).'* ar1 * b;
        
        %scaling alpha pass
        format long
        s(1,k)= alpha(:,k).'*one;
        format long
        salphainter(:,k)= alpha(:,k)./s(1,k);
        
        if salphainter(1,k)<1e-135 || isnan(salphainter(1,k))
            salphainter(1,k)=0;
        elseif salphainter(2,k)<1e-135|| isnan(salphainter(2,k))
            salphainter(2,k)=0;
        end
        salpha(:,k)= salphainter(:,k);
    end
    
    %% Beta pass
    format long
    s_b = 1./ s;
    format long
    sbeta(:,T)=[s_b(1,T);s_b(1,T)];
    
    for j = T-1:-1:1
        format long
        f_1 =((mu(1,1)*o(j+1,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(o(j+1,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
        f_2 =((mu(1,2)*o(j+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(j+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
        b(1,1)= f_1;
        b(2,2)= f_2;
        format long
        beta(:,j)= ar1*b*sbeta(:,j+1);
        %beta scaling
        format long
        sbeta(:,j)= beta(:,j).*s_b(1,j);
    end
    format long
    sbeta(1,T) = 1./s(1,T);
    
    %% gamma
    
    for k = 1:(T-1)
        format long
        f_1 =((mu(1,1)*o(k+1,1)^(mu(1,1)*p(1,1)-1))/(sigma(1,1)^(mu(1,1)*p(1,1))*(gamma(p(1,1))*gamma(q(1,1))/gamma(p(1,1)+q(1,1)))*(1+(o(k+1,1)/sigma(1,1))^mu(1,1))^(p(1,1)+q(1,1))));
        f_2 =((mu(1,2)*o(k+1,1)^(mu(1,2)*p(1,2)-1))/(sigma(1,2)^(mu(1,2)*p(1,2))*(gamma(p(1,2))*gamma(q(1,2))/gamma(p(1,2)+q(1,2)))*(1+(o(k+1,1)/sigma(1,2))^mu(1,2))^(p(1,2)+q(1,2))));
        
        b(1,1)= f_1;
        b(2,2)= f_2;
        
        format long
        sa_1 = [salpha(1,k),0;0,salpha(1,k)];
        
        
        matb = vec2mat(sbeta(:,k+1),4);
        matb_1 = reshape(matb,2,2);
        matb_1(2,2)=matb_1(2,1);
        matb_1(2,1)= 0;
        
        tra = vec2mat(ar1(1,:),4);
        tra_1 = reshape (tra,2,2);
        tra_1(2,2) = tra_1(2,1);
        tra_1(2,1)=0;
        
        format long
        digam = sa_1*matb_1*tra_1*b*one;
        
        
        gammma(1,k)=(one.'*digam);
        
        
        sa_2 = [salpha(2,k),0;0,salpha(2,k)];
        
        tra2 = vec2mat(ar1(2,:),4);
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
    % gamma
    gammma(:,T)= alpha(:,T)./ s(1,T);    
    %% calcule P(O/lambda)
    Pb(1,L) = -((log10(s))*ones(T,1))
    
    %% restimation
    
    %pi
    pi1= gammma(:,1).';
    
    % A
    
    for e = 1:2
        ar1(1,e)= (digamma(e,:)*ones(T-1,1)) / (gammma(1,(1:T-1))*ones(T-1,1));
        ar1(2,e)= (digamma(e+2,:)*ones(T-1,1))/(gammma(2,(1:T-1))*ones(T-1,1));
    end
    
    
    % B
    f1=@(x)-sum(sum(gammma(1,:),2)*(log(x(1))-log(x(2))-log(gamma(x(3)))-log(gamma(x(4)))+log(gamma(x(3)+x(4))))+gammma(1,:)*(x(1)*x(3)-1)*log(o/x(2))-gammma(1,:)*(x(3)+x(4))*log(1+(o/x(2)).^x(1)));
    lb = zeros(1,4);
    ub = +Inf;
    xfinal1 = fmincon(f1,[muini(1,1),sigmaini(1,1),pini(1,1),qini(1,1)],[],[],[],[],lb,ub)
    
    
    
    f2=@(x)-sum(sum(gammma(2,:),2)*(log(x(1))-log(x(2))-log(gamma(x(3)))-log(gamma(x(4)))+log(gamma(x(3)+x(4))))+gammma(2,:)*(x(1)*x(3)-1)*log(o/x(2))-gammma(2,:)*(x(3)+x(4))*log(1+(o/x(2)).^x(1)));
    lb = zeros(1,4);
    ub = +Inf;
    xfinal2 = fmincon(f2,[muini(1,2),sigmaini(1,2),pini(1,2),qini(1,2)],[],[],[],[],lb,ub)
    
    mu = [xfinal1(1,1),xfinal2(1,1)];
    sigma = [(xfinal1(1,2)),(xfinal2(1,2))];
    p = [xfinal1(1,3),xfinal2(1,3)];
    q = [xfinal1(1,4),xfinal2(1,4)];
    
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
iteration =  ones(1,length(Pb));
for l = 1 : length(Pb)-1
    iteration(1,l+1)=iteration(1,l)+1;
end    
proba = plot(iteration,Pb(1,:));
drawnow
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('P(O/lambda)')
xlabel('Iterations')
title('convergence des proba','FontSize',8)
hold off

%% graf de verification 

%densité gb2
a1= xfinal1(1,1);
b1 =xfinal1(1,2);
p1= xfinal1(1,3);
q1= xfinal1(1,4);
a2= xfinal2(1,1);
b2 =xfinal2(1,2);
p2= xfinal2(1,3);
q2= xfinal2(1,4);

funct1 =@(z)(taille(1,1)*((a1*z.^(a1*p1-1))./(b1^(a1*p1)*(gamma(p1)*gamma(q1)/gamma(p1+q1))*(1+(z/b1).^a1).^(p1+q1)))...
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

% graf des densites
figure(figure_2)
subplot(2,2,1)
hold on
ezplot(@(x)(taille(1,1)*((a1*x.^(a1*p1-1))./(b1^(a1*p1)*(gamma(p1)*gamma(q1)/gamma(p1+q1))*(1+(x/b1).^a1).^(p1+q1)))...
            +taille(2,1)*((a2*x.^(a2*p2-1))./(b2^(a2*p2)*(gamma(p2)*gamma(q2)/gamma(p2+q2))*(1+(x/b2).^a2).^(p2+q2)))),[0,1000,0,0.01]);
title('Fonction de densité', 'FontSize', 8)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('density')
xlabel('observation en million')
hold off


%% influence des PDF Les deux CDF sont comparés   

GB2_1=zeros(T,1);
for t=1:T
    GB2_1(t,1) = GB2CDF(o(t,1),mu(1,1),sigma(1,1),p(1,1),q(1,1));
end

GB2_2 = zeros(T,1);
for t=1:T
    GB2_2(t,1) = GB2CDF(o(t,1),mu(1,2),sigma(1,2),p(1,2),q(1,2));
end

testpdf=zeros(T,1);
for i=1:T
    testpdf(i,1) = taille(2,1)*GB2_2(i,1)+taille(1,1)*GB2_1(i,1);
end

figure(figure_2)
subplot(2,2,2)
hold on
scatter(o,testpdf,10,'o','filled')
scatter(o,GB2_1,'.')
scatter(o,GB2_2,5,'d')
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
parampdf =@(vari) (taille(1,1)*(GB2CDF(vari,mu(1,1),sigma(1,1),p(1,1),q(1,1)))...
    + taille(2,1)*GB2CDF(vari,mu(1,2),sigma(1,2),p(1,2),q(1,2)));

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

AIC = -2*(-Pb(1,L))+22;
BIC = -2*(-Pb(1,L))+11*log(T);

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
    EmQ(i,1)= (1/taille(1,1))*(GB2qtl(Quantiltest(i,1),mu(1,1),sigma(1,1),p(1,1),q(1,1)))...
      +(1/taille(2,1))*(GB2qtl(Quantiltest(i,1),mu(1,2),sigma(1,2),p(1,2),q(1,2)));
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

%% creation de la figure des gammabis
figure_3 = figure
t =(1:1:length(o));
hold on
plot(t,gammabis2(:,2),'--k')
hold on
drawnow; 
legend('régime 2')
title('régime', 'FontSize',8)
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
A1 = [mu(1,1);mu(1,2);A(1,1);A(1,2)];
B1 = [sigma(1,1);sigma(1,2);B(1,1);B(1,2)];
P1 = [p(1,1);p(1,2);P(1,1);P(1,2)];
Q1 = [q(1,1);q(1,2);Q(1,1);Q(1,2)];
table_parametres = table(A1,B1,P1,Q1,'RowNames',parametre)

figure_5 = figure;
format short
t = uitable(figure_5,'Data',{'Loi1final',mu(1,1),sigma(1,1),p(1,1),q(1,1)...
    ;'loi2final',mu(1,2),sigma(1,2),p(1,2),q(1,2);'loi1ini',A(1,1),...
    B(1,1),P(1,1),Q(1,1);'loi2ini',A(1,2),...
    B(1,2),P(1,2),Q(1,2)},'Position',[85 250 410 102]);
t.ColumnName = {'loi','A','B','P','Q'};
t.ColumnEditable = true;
ax = gca;
set(ax, 'Visible','off')
hold off

%% creation des tablo crise non crise
pbnocriseL1=(sum(gammabis2(1:27,1))+sum(gammabis2(36:107,1))...
    +sum(gammabis2(127:150,1)))/(150-27);

pbnocriseL2=(sum(gammabis2(1:27,2))+sum(gammabis2(36:107,2))...
    +sum(gammabis2(127:150,2)))/(150-27);

pbcriseL1=(sum(gammabis2(28:35,1))+sum(gammabis2(108:126,1)))/27;

pbcriseL2=(sum(gammabis2(28:35,2))+sum(gammabis2(108:126,2)))/27;

format short
tete = {'crise';'non-crise'};
Loi_1 = [pbcriseL1;pbnocriseL1];
Loi_2 = [pbcriseL2;pbnocriseL2];
table_crise = table(Loi_1,Loi_2,'RowNames',tete)


figure(figure_5);
format short
cris = uitable(figure_5,'Data',{'Loi_1',pbcriseL1,pbnocriseL1...
    ;'Loi_2',pbcriseL2,pbnocriseL2},'Position',[85 120 259 75]);
cris.ColumnName = {'loi','Crise','Non Crise'};
cris.ColumnEditable = true;
ax = gca;
set(ax, 'Visible','off')
hold off

print('figure_5','-djpeg')
print ('figure_5','-dpdf','-r600')
mydoc;

