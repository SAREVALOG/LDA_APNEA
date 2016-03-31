% Carga candidatos a clasificar de una base de datos en EXCEL
% grafica espacios de dispersión 
% usa LDA para clasificar 
%
%      created on 2015 by
%      Santiago Arevalo (National University of Entre Rios)
%      (email: sarevalog@correo.udistrital.edu.co)
% -------------------------------------------------------------------------

clc; clear all;

PATH= '/home/sarevalog/Dropbox/MAESTRIA/THESIS/CODIGOS/MATLAB/RESPIRATORY';

MATRIZ = xlsread('descriptores.xlsx','final');

M= MATRIZ(:,1);
MIN= MATRIZ(:,2);
LESS= MATRIZ(:,3);
enerVLF=MATRIZ(:,4);
RMSSD=MATRIZ(:,5);
std_RR=MATRIZ(:,6);
mean_RR=MATRIZ(:,7);
apnea=MATRIZ(:,9);

fprintf(1,'\\n$> LOAD FINISH SUCCESSFUL \n');

% %% --------- Graficar espacio de caracteristicas de dos 
% %---------- dimensiones (M , MIN) 
% 
n=length(M);
y1=apnea(1:n);
mean_spo2=M(1:n);
min_spo2=MIN(1:n);
less_spo2=LESS(1:n);
enerVLF=enerVLF(1:n);
mean_RR=mean_RR(1:n);
std_RR=std_RR(1:n);
RMSSD = RMSSD(1:n);

% ---------- dimensiones (MIN , M) 
X1=[min_spo2,mean_spo2];
figure(1), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea','Location','northwest');
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('Minimo valor de SpO2 (%)');



% ---------- dimensiones (MIN , LESS) 
% 
 X1=[min_spo2,less_spo2];
 figure(2),  clf, box on, hold on, grid on;
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('Numero de valores menores al 92 en 1 minuto');
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea','Location','northwest');

% plot(X1(1,y1==1),X1(2,y1==1),'r+',X1(1,y1==8),X1(2,y1==8),'bo')

%---------- dimensiones (MIN , STDRR) 
 
 X1=[min_spo2,std_RR];
 figure(3),  clf, box on, hold on, grid on;
 xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('Numero de valores menores al 92 en 1 minuto');
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
% plot(X1(1,y1==1),X1(2,y1==1),'r+',X1(1,y1==8),X1(2,y1==8),'bo')


%-------------------------(minspo2, EnerVLF)

X1=[min_spo2,enerVLF];
 
figure(4), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Desviación típica de los picos RR (ms)'); 
 ylabel('FMB (ms^2)');

%-------------------------(minspo2, mean(RR))
X1=[min_spo2,mean_RR];
figure(5), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Desviación típica de los picos RR (ms)'); 
 ylabel('Valor medio de los picos RR (ms)');

% -------------------------(minspo2, RMSSD)
X1=[min_spo2, RMSSD];
figure(6), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('FMB (ms^2)'); 
 ylabel('Valor medio de los picos RR (ms)');


%-------------------------(meam_spo2,less)
X1=[mean_spo2, less_spo2];
figure(7), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('Desviación típica de los picos RR (ms)');

%-------------------------(mean_spo2,std_RR)
X1=[mean_spo2, std_RR];
figure(8), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('FMB (ms^2)');

%-------------------------(mean_spo2,enerVLF)
X1=[mean_spo2, enerVLF];
figure(9), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('Valor medio de los picos RR (ms)');

%-------------------------(mean_spo2,meanRR)
X1=[mean_spo2, mean_RR];
figure(10), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('Desviación típica de los picos RR (ms)');

%-------------------------(mean_spo2,RMSSD)
X1=[mean_spo2, RMSSD];
figure(11), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('FMB (ms^2)');

%-------------------------(less_spo2,std_RR)
X1=[less_spo2, std_RR];
figure(12), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('Valor medio de los picos RR (ms)');
 
%-------------------------(less_spo2,ener_VLF)
X1=[less_spo2, enerVLF];
figure(13), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Numero de valores menores al 92 en 1 minuto'); 
 ylabel('Desviación típica de los picos RR (ms)');
 
%-------------------------(less_spo2,meanRR)
X1=[less_spo2, mean_RR];
figure(14), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Numero de valores menores al 92 en 1 minuto'); 
 ylabel('FMB (ms^2)');
 
%-------------------------(less_spo2,rmssd)
X1=[less_spo2, RMSSD];
figure(15), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('Numero de valores menores al 92 en 1 minuto'); 
 ylabel('Valor medio de los picos RR (ms)');
 
%--------------------------(stdRR,enerVLF)
X1=[std_RR, enerVLF];
figure(16), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('desviación típica RR (ms)');

%--------------------------(stdRR,meanRR)
X1=[std_RR, mean_RR];
figure(17), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('FMB (ms^2)');
 
%--------------------------(stdRR,RMSSD)
X1=[std_RR, RMSSD];
figure(18), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('promedio picos RR (ms)');
 

%--------------------------(enerVLF, meanRR)
X1=[enerVLF, mean_RR];
figure(19), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('Numero de valores menores al 92 en 1 minuto');


%--------------------------(enerVLF,RMSSD)
X1=[enerVLF, RMSSD];
figure(20), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('promedio SpO_2 (%)');

 
%--------------------------(mean_RR,RMSSD)
X1=[mean_RR, RMSSD];
figure(21), clf, box on, hold on, grid on;
gscatter(X1(:,1),X1(:,2),y1,'br','xo');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('minimo valor SpO2 (%)');
 
 
fprintf(1,'\\n$> PLOT FEATURES FINISHED \n');

%% -------------VERIFICACIÓN DE NORMALIDAD -----
h_minspo2=lillietest(min_spo2');

[f,x_values] = ecdf(min_spo2');
J = plot(x_values,f);
hold on;
K = plot(x_values,normcdf(x_values),'r--');
set(J,'LineWidth',2);
set(K,'LineWidth',2);
legend([J K],'min_spo2','Standard Normal CDF','Location','SE');

h_meanspo2=lillietest(mean_spo2');
h_lessspo2=lillietest(less_spo2');
h_stdRR=lillietest(std_RR');
h_enerVLF=lillietest(enerVLF');
h_meanRR=lillietest(mean_RR');
h_RMSSD=lillietest(RMSSD');


fprintf(1,'\\n$> KOLMOGOROV-SMIRNOV TEST FINISHED \n');


% -------------DISCRIMINANTE LINEAL ----------------------------
%------------- evaluar el comportamiento de los 
%------------- candidatos
X1=[min_spo2];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr_a = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM_a,grpOrder_a] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres_a] = confusionmatStats(y1,ldaClass); % Evaluacion

X1=[mean_spo2];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr_b = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM_b,grpOrder_b] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres_b] = confusionmatStats(y1,ldaClass); % Evaluacion

X1=[less_spo2];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr_c = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM_c,grpOrder_c] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres_c] = confusionmatStats(y1,ldaClass); % Evaluacion

X1=[enerVLF];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr_d = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM_d,grpOrder_d] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres_d] = confusionmatStats(y1,ldaClass); % Evaluacion

X1=[RMSSD];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr_e = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM_e,grpOrder_e] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres_e] = confusionmatStats(y1,ldaClass); % Evaluacion

X1=[std_RR];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr_f = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM_f,grpOrder_f] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres_f] = confusionmatStats(y1,ldaClass); % Evaluacion

X1=[mean_RR];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr_g = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM_g,grpOrder_g] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres_g] = confusionmatStats(y1,ldaClass); % Evaluacion

%--------------------------espacios dos dimensiones
% ---------- dimensiones (MIN , M) 
X1=[min_spo2,mean_spo2];

lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr1 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM1,grpOrder1] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres1] = confusionmatStats(y1,ldaClass); % Evaluacion

figure (1),
%Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 110 0 90]);
h2.Color = 'r';
h2.LineWidth = 2;
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('Minimo valor de SpO2 (%)');
title('{\bf Clasificación con LDA}');

% ---------- dimensiones (MIN , LESS) 
% 
 X1=[min_spo2,less_spo2];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr2 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM2,grpOrder2] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres2] = confusionmatStats(y1,ldaClass);

figure (2),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 110 0 7000]);
h2.Color = 'r';
h2.LineWidth = 2;
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('Numero de valores menores al 92 en 1 minuto');
 title('{\bf Clasificación con LDA}')

%---------- dimensiones (MIN , STDRR) 
 
 X1=[min_spo2,std_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr3 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM3,grpOrder3] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres3] = confusionmatStats(y1,ldaClass);

figure (3),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 100 0 7000]);
h2.Color = 'k';
h2.LineWidth = 2;
xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('Numero de valores menores al 92 en 1 minuto');
 title('{\bf Clasificación con LDA}');
 

%-------------------------(minspo2, EnerVLF)

X1=[min_spo2,enerVLF];
 
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr4 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM4,grpOrder4] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres4] = confusionmatStats(y1,ldaClass);

figure (4),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-20 2000 0 7000]);
h2.Color = 'k';
h2.LineWidth = 2;
 xlabel('Desviación típica de los picos RR (ms)'); 
 ylabel('FMB (ms²)');
 title('{\bf Clasificación con LDA}');

%-------------------------(minspo2, mean(RR))
X1=[min_spo2,mean_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr5 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM5,grpOrder5] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres5] = confusionmatStats(y1,ldaClass);

figure (5),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 2000 0 4000]);
h2.Color = 'k';
h2.LineWidth = 2;
 xlabel('Desviación típica de los picos RR (ms)'); 
 ylabel('Valor medio de los picos RR (ms)');
 title('{\bf Clasificación con LDA}');
 

% -------------------------(minspo2, RMSSD)
X1=[min_spo2, RMSSD];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr6 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM6,grpOrder6] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres6] = confusionmatStats(y1,ldaClass);

figure (6),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-2000 7000 0 7000]);
h2.Color = 'k';
h2.LineWidth = 2;
 xlabel('FMB (ms²)'); 
 ylabel('Valor medio de los picos RR (ms)');
 title('{\bf Clasificación con LDA}');

%-------------------------(mean_spo2,less)
X1=[mean_spo2, less_spo2];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr7 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM7,grpOrder7] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres7] = confusionmatStats(y1,ldaClass);

figure (7),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 100 0 2000]);
h2.Color = 'k';
h2.LineWidth = 2;
xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('Desviación típica de los picos RR (ms)');
 title('{\bf Clasificación con LDA}');
 
%-------------------------(mean_spo2,std_RR)
X1=[mean_spo2, std_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr8 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM8,grpOrder8] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres8] = confusionmatStats(y1,ldaClass);

figure (8),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 100 0 70000]);
h2.Color = 'k';
h2.LineWidth = 2;
 xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('FMB (ms²)');
 title('{\bf Clasificación con LDA}');

%-------------------------(mean_spo2,enerVLF)
X1=[mean_spo2, enerVLF];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr9 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM9,grpOrder9] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres9] = confusionmatStats(y1,ldaClass);

figure (9),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 100 0 3500]);
h2.Color = 'k';
h2.LineWidth = 2;
 xlabel('Minimo valor de SpO2 (%)'); 
 ylabel('Valor medio de los picos RR (ms)');
 title('{\bf Clasificación con LDA}');

%-------------------------(mean_spo2,meanRR)
X1=[mean_spo2, mean_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr10 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM10,grpOrder10] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres10] = confusionmatStats(y1,ldaClass);

figure (10),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 100 0 2000]);
h2.Color = 'k';
h2.LineWidth = 2;
xlabel('Valor medio de SpO2 (%)'); 
ylabel('Desviación típica de los picos RR (ms)');
title('{\bf Clasificación con LDA}');
%-------------------------(mean_spo2,RMSSD)
X1=[mean_spo2, RMSSD];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr11 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM11,grpOrder11] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres11] = confusionmatStats(y1,ldaClass);

figure (11),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 100 0 70000]);
h2.Color = 'k';
h2.LineWidth = 2;
xlabel('Valor medio de SpO2 (%)'); 
 ylabel('FMB (ms²)');
title('{\bf Clasificación con LDA}');

%-------------------------(less_spo2,std_RR)
X1=[less_spo2, std_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr12 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM12,grpOrder12] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres12] = confusionmatStats(y1,ldaClass);

figure (12),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 100 0 3500]);
h2.Color = 'k';
h2.LineWidth = 2;
 xlabel('Valor medio de SpO2 (%)'); 
 ylabel('Valor medio de los picos RR (ms)');
title('{\bf Clasificación con LDA}');

%-------------------------(less_spo2,ener_VLF)
X1=[less_spo2, enerVLF];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr13 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM13,grpOrder13] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres13] = confusionmatStats(y1,ldaClass);

figure (13),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 6000 0 2000]);
h2.Color = 'k';
h2.LineWidth = 2;
xlabel('Numero de valores menores al 92 en 1 minuto'); 
 ylabel('Desviación típica de los picos RR (ms)');
 title('{\bf Clasificación con LDA}');
 
%-------------------------(less_spo2,meanRR)
X1=[less_spo2, mean_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr14 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM14,grpOrder14] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres14] = confusionmatStats(y1,ldaClass);

figure (14),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 7000 0 70000]);
h2.Color = 'k';
h2.LineWidth = 2;
xlabel('Numero de valores menores al 92 en 1 minuto'); 
 ylabel('FMB (ms²)');
  title('{\bf Clasificación con LDA}');

 %-------------------------(less_spo2,rmssd)
X1=[less_spo2, RMSSD];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr15 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM15,grpOrder15] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres15] = confusionmatStats(y1,ldaClass);

figure (15),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 2000 0 3500]);
h2.Color = 'k';
h2.LineWidth = 2;
xlabel('Numero de valores menores al 92 en 1 minuto'); 
 ylabel('Valor medio de los picos RR (ms)');
   title('{\bf Clasificación con LDA}');
   
%--------------------------(stdRR,enerVLF)
X1=[std_RR, enerVLF];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr16 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM16,grpOrder16] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres16] = confusionmatStats(y1,ldaClass);

figure (16),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 1500 0 2500]);
h2.Color = 'k';
h2.LineWidth = 2;
title('{\bf Clasificación con LDA}');
xlabel('RMSSD (ms)'); 
 ylabel('Desviación típica RR (ms)');
   title('{\bf Clasificación con LDA}');

%--------------------------(stdRR,meanRR)
X1=[std_RR, mean_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr17 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM17,grpOrder17] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres17] = confusionmatStats(y1,ldaClass);

figure (17),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 1500 0 3500]);
h2.Color = 'k';
h2.LineWidth = 2;
title('{\bf Clasificación con LDA}');
xlabel('RMSSD (ms)'); 
 ylabel('FMB (ms²)');
   title('{\bf Clasificación con LDA}');

%--------------------------(stdRR,RMSSD)
X1=[std_RR, RMSSD];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr18 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM18,grpOrder18] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres18] = confusionmatStats(y1,ldaClass);

figure (18),
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[-10 1500 0 3500]);
h2.Color = 'k';
h2.LineWidth = 2;
title('{\bf Clasificación con LDA}');
xlabel('RMSSD (ms)'); 
 ylabel('Promedio intervalo RR (ms)');
   title('{\bf Clasificación con LDA}');

%--------------------------(enerVLF, meanRR)
X1=[enerVLF, mean_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr19 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM19,grpOrder19] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres19] = confusionmatStats(y1,ldaClass);

figure(19), 
% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 6 0 6000]);
h2.Color = 'k';
h2.LineWidth = 2;
   title('{\bf Clasificación con LDA}');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('Numero de valores menores al 92 en 1 minuto');


%--------------------------(enerVLF,RMSSD)
X1=[enerVLF, RMSSD];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr20 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM20,grpOrder20] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres20] = confusionmatStats(y1,ldaClass);

figure(20), 
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 6 0 100]);
h2.Color = 'k';
h2.LineWidth = 2;
   title('{\bf Clasificación con LDA}');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('promedio SpO_2 (%)');

 
%--------------------------(mean_RR,RMSSD)
X1=[mean_RR, RMSSD];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErr21 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCM21,grpOrder21] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatres21] = confusionmatStats(y1,ldaClass);

figure(21), 
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2=ezplot(f,[0 6 0 100]);
h2.Color = 'k';
h2.LineWidth = 2;
   title('{\bf Clasificación con LDA}');
legend('normal','apnea');
 xlabel('RMSSD (ms)'); 
 ylabel('minimo valor SpO2 (%)');
 
%% ---------------------EVALUACION MEJOR DESCRIPTOR -------------

X1=[min_spo2, less_spo2, mean_RR];
lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
K_x1= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L_x1= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErrfinal1 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCMfinal1,grpOrderfinal1] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatresfinal1] = confusionmatStats(y1,ldaClass);

cvmodel1 = crossval(lda,'kfold',10);
cverror1 = kfoldLoss(cvmodel1);

X2=[min_spo2, less_spo2, enerVLF];
lda = fitcdiscr(X2,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErrfinal2 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCMfinal2,grpOrderfinal2] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatresfinal2] = confusionmatStats(y1,ldaClass);

cvmodel2 = crossval(lda,'kfold',10);
cverror2 = kfoldLoss(cvmodel2);

X3=[min_spo2, less_spo2, enerVLF, mean_RR];
lda = fitcdiscr(X3,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErrfinal3 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCMfinal3,grpOrderfinal3] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatresfinal3] = confusionmatStats(y1,ldaClass);

cvmodel3 = crossval(lda,'kfold',10);
cverror3 = kfoldLoss(cvmodel3);

X4=[min_spo2, mean_spo2, less_spo2, std_RR, enerVLF, mean_RR, RMSSD];
lda = fitcdiscr(X4,y1);             % Se crea clasificador lineal
K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
ldaResubErrfinal4 = resubLoss(lda);      % Resubstitution error 
ldaClass = resubPredict(lda);
[ldaResubCMfinal4,grpOrderfinal4] = confusionmat(y1,ldaClass); % confusion matrix
[confusionmatresfinal4] = confusionmatStats(y1,ldaClass);

cvmodel4 = crossval(lda,'kfold',10);
cverror4 = kfoldLoss(cvmodel4);

