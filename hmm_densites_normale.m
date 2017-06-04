clear ''all'';
close all;
clc;
%%
%%%%% Paramètres initiaux
% 'a' est la matrice de transition initiale 
% 'moyenne' est le vecteur de moyennes initiales des deux densitésnormales 
% 'etype' est le vecteur d'écart-types des deux lois normales initiales
% 'pi' est la probabilité initiale régime 1 ou régime 2 
% N est le nombre de régimes utilisés
% maxit est le nombre d'itérations souhaité pour le modèle

a = [0.50,0.50;0.50,0.50];

moyenne=[0.1,50];
etype=[500,500];
pi1= [0.4,0.6];
maxit = 120;

mu=[moyenne(1,1),moyenne(1,2)];
sigma = [etype(1,1),etype(1,2)];

% Importation des données
o = xlsread('tablo time series.xlsx','TIME_SERIES','D2:D151');

%%% création des matrices des poids
T = length(o);
gamma = ones(N,T);
digamma= ones(2*N,T-1);
one = [1;1];
N=2;

%% definition du temps 
tps(1,1)=1;
for t = 1:T-1
    tps(t+1,1)=tps(t,1)+1;
end
%% Algorithme Baum_Welch

for L= 1:maxit
    %% alpha pass
    f_1 = normpdf(o(1,1),mu(1,1),sigma(1,1));
    f_2 = normpdf(o(1,1),mu(1,2),sigma(1,2));

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

    for k=2:T
        format long
        f_1 = normpdf(o(k,1),mu(1,1),sigma(1,1));
        f_2 = normpdf(o(k,1),mu(1,2),sigma(1,2));
        b(1,1)= f_1;
        b(2,2)= f_2; 
        alpha(:,k) = salpha(:,k-1).'* a * b;
 
    %scaling alpha pass
        format long  
        s(1,k)= alpha(:,k).'*one;
        format long
        salpha(:,k)= alpha(:,k)./s(1,k);
    end

    %% Beta pass
    format long
    s_b = 1./ s;
    format long
    sbeta(:,T)=[s_b(1,T);s_b(1,T)];

    for j = T-1:-1:1
        format long
        f_1 = normpdf(o(j+1,1),mu(1,1),sigma(1,1));
        f_2 = normpdf(o(j+1,1),mu(1,2),sigma(1,2));
        b(1,1)= f_1;
        b(2,2)= f_2;
        format long
        beta(:,j)= a*b*sbeta(:,j+1);
    %beta scaling
        format long
        sbeta(:,j)= beta(:,j).*s_b(1,j);
    end
    format long
    sbeta(:,T) = 1./s(1,T);

    %% gamma

    for k = 1:(T-1)
        format long
        f_1 = normpdf(o(k+1,1),mu(1,1),sigma(1,1));
        f_2 = normpdf(o(k+1,1),mu(1,2),sigma(1,2));
    
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


        gamma(1,k)=(one.'*digam);
    
 
        sa_2 = [salpha(2,k),0;0,salpha(2,k)];

        tra2 = vec2mat(a(2,:),4);
        tra_2 = reshape (tra2,2,2);
        tra_2(2,2) = tra_2(2,1);
        tra_2(2,1)=0;

        format long
        digam2 =  sa_2*matb_1*tra_2*b*one;
        gamma(2,k)= (one.'*digam2);
       
        % reshaping digamma
        digam_1 = vec2mat(digam,4,[digam2(1,1),digam2(2,1)]);
        digamma(:,k) = reshape(digam_1,4,1);

    end
    % gamma_T
    gamma(:,T)= alpha(:,T)./ s(1,T);
    %% calcule P(O/lambda)
    P(1,L) = -((log10(s))*ones(T,1))
    
    %% restimation
            pi1= gamma(:,1)';
            
            % A
            
            for e = 1:2
                a(1,e)= (digamma(e,:)*ones(T-1,1)) / (gamma(1,(1:T-1))*ones(T-1,1));
                a(2,e)= (digamma(e+2,:)*ones(T-1,1))/(gamma(2,(1:T-1))*ones(T-1,1));
            end
            
            % B
            %variance
            sigma_1=[0,0];
            
            for t = 1:T
                sigma_1(1,1)= (((o(t,1)-mu(1,1))^2*gamma(1,t)))+sigma_1(1,1);
                sigma_1(1,2)= (((o(t,1)-mu(1,2))^2*gamma(2,t)))+ sigma_1(1,2);
            end
            
            sigma(1,1)=sqrt(sigma_1(1,1)/(gamma(1,:)*ones(T,1)));
            sigma(1,2)=sqrt(sigma_1(1,2)/(gamma(2,:)*ones(T,1)));
   
            
            %moyenne
            mu(1,1)= (gamma(1,:)*o)/(gamma(1,:)*ones(T,1));
            mu(1,2)= (gamma(2,:)*o)/(gamma(2,:)*ones(T,1));
           
     %% graphique
     if L~=maxit && L~=1
         if mod(L,10)==0 
             figure;
             hold on
             art = 25;
             cc = zeros(T,3);
             gammabis2 = gamma';
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
         gammabis2 = gamma';
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
         set(gca, 'FontName', 'Arial')
         set(gca, 'FontSize', 8)
         ylabel('Log(data)')
         xlabel('Temps')
         set(gca,'XTick',[]);
         caxis([0 1])
         drawnow
         title('Observations: dernière itération','FontSize',8)
         hold off
     elseif L==1
         figure_2=figure ;
         subplot(2,2,3)
         hold on
         art = 25;
         cc = zeros(T,3);
         gammabis2 = gamma';
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
%%
%%%%%%%%%%%%%%%%%%%%%%%% Rapport de fin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calcule de la taille 
taille = [(gamma(1,:)*ones(T,1)/length(o)); ((gamma(2,:)*ones(T,1))/length(o))];

%% graf de la convergences des probabilités
figure_4 = figure;
hold on
iteration =  ones(1,maxit);
for l = 1 : maxit-1
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


% graphique des densités estimées
figure;
hold on
sigma_graf = cat(3, sigma(1,1), sigma(1,2));
gm = gmdistribution(mu', sigma_graf, taille);
gmPDF = @(x)pdf(gm,[x]);
ezplot(gmPDF,[0 1000]);
drawnow
hold off

% graphique densités histogramme
for i = 1:T
    found(i,1)=taille(1,1)*normpdf(o(i,1),mu(1,1),sigma(1,1))+...
        taille(2,1)*normpdf(o(i,1),mu(1,2),sigma(1,2));
end

figure_1 = figure;
hold on 
subplot(2,2,1)
histogram((o),150,'Normalization','pdf','FaceColor','b')
hold on
h=scatter((o),found,5,'filled','m');
axis([-inf,inf,0,0.01]);
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('density')
xlabel('observation en million')
hold off
title('Donnée et densité', 'FontSize', 8)

% graphique mixte des densites
figure(figure_2)
subplot(2,2,1)
fplot(@(x)taille(1,1)*normpdf(x,mu(1,1),sigma(1,1))...
    +taille(2,1)*normpdf(x,mu(1,2),sigma(1,2)),[0,max(o)]);
title('Fonction de densité', 'FontSize', 8)
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
ylabel('density')
xlabel('observation en million')
hold off


% influence des PDF Les deux CDF sont comparées   

GB2_1=zeros(T,1);
for t=1:T
    GB2_1(t,1) = normcdf(o(t,1),mu(1,1),sigma(1,1));
end

GB2_2 = zeros(T,1);
for t=1:T
    GB2_2(t,1) = normcdf(o(t,1),mu(1,2),sigma(1,2));
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


% comparaison PDF et test  
parampdf =@(vari) (taille(1,1)*(normcdf(vari,mu(1,1),sigma(1,1)))...
    + taille(2,1)*normcdf(vari,mu(1,2),sigma(1,2)));

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test de bon fit du modèle avec pseudo residus 

unipseudo=zeros(T,1);
for i = 1:T
    unipseudo(i,1)= parampdf(o(i,1));   
end
npseudores = zeros(T,1);
for i=1:T
    npseudores(i,1) = norminv(unipseudo(i,1),0,1);
end

% normal QQ plot 
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

AIC = -2*log(P(1,L))+2*7;
BIC = -2*log(P(1,L))+7*log(T);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creation des QQ plots pour les donnees 
Quantiltest = [0.01:0.05:0.99]';
EmQ=zeros(length(Quantiltest),1);

for i=1:length(Quantiltest)
    EmQ(i,1)= (1/taille(1,1))*(norminv(Quantiltest(i,1),mu(1,1),sigma(1,1)))...
      +(1/taille(2,1))*(norminv(Quantiltest(i,1),mu(1,2),sigma(1,2)));
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
plot(t,gamma(1,:),'k')
%hold on 
%plot(t,gamma(2,:),':r')
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
A = [mu(1,1);mu(1,2);moyenne(1,1);moyenne(1,2)];
B = [sigma(1,1);sigma(1,2);etype(1,1);etype(1,2)];
table_parametres = table(A,B,'RowNames',parametre)

figure_5 = figure;
format short
t = uitable(figure_5,'Data',{'Loi1final',mu(1,1),sigma(1,1)...
    ;'loi2final',mu(1,2),sigma(1,2);'loi1ini',moyenne(1,1),...
    etype(1,1);'loi2ini',moyenne(1,2),...
    etype(1,2)},'Position',[85 250 410 102]);
t.ColumnName = {'loi','Mu','Sigma'};
t.ColumnEditable = true;
ax = gca;
set(ax, 'Visible','off')
hold off

print('figure_5','-djpeg')
print ('figure_5','-dpdf','-r600')

% Rapport sur word
mydoc;
