% Este programa lee los datos de ECG los cuales estan guardados con el 
% formato 212 y datos de oximetria y respiración guardados en formato 16
% de la base de datos APNEA-ECG disponibles en physionet.org
% Las anotaciones han sido preparadas basadas en la revisión visual de 
% un experto revisando simultaniamente las señales de respiracion y 
% saturación de oxigeno.
% Cuando estos archivos fueron publicados, la siguiente (incorrecta)
% descripción de apnean fue usada:
% En estos archivos, un "8" indicua que ocurre un evento apenico durante el
% siguiente intervalo de un minuto, y un "1" indicua que no existe apnea
% durante el siguiente intervalo de un minuto. 
%
%
%      created on 2015 by
%      Santiago Arevalo (National University of Entre Rios)
%      (email: sarevalog@correo.udistrital.edu.co),
%
%      algorithm is based on a program written by
%      Klaus Rheinberger (University of Innsbruck)
%      (email: klaus.rheinberger@uibk.ac.at)
%
% -------------------------------------------------------------------------
clc; clear all;

PATH= '/home/sarevalog/Dropbox/MAESTRIA/THESIS/CODIGOS/MATLAB/RESPIRATORY';
sig='c01';

DATAFILEECG=[sig '.dat'];         % data-file
DATAFILERES=[sig 'r' '.dat'];         % data-file
ATRFILE = [sig 'r' '.apn'];       % apnea annotations
HEADERFILE= [sig 'er' '.hea'];      % header-file in text format

%% ------ CARGA DATOS BINARIOS ----------------------------------------------
% ------ CARGA ENCABEZADO DE LOS DATOS ------------------------------------

fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE);
signalh= fullfile(PATH, HEADERFILE);
fid1=fopen(signalh,'r');
z= fgetl(fid1);
A= sscanf(z, '%*s %d %d %d',[1,3]);
nosig= A(1);  % number of signals
sfreq=A(2);   % sample rate of data
samples = A(3);
clear A;
for k=1:nosig
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d %d %d %*s',[1,5]);
    dformat(k)= A(1);           % format; 
    gain(k)= A(2);              % number of integers per mV
    bitres(k)= A(3);            % bitresolution
    zerovalue(k)= A(4);         % integer value of ECG zero point
    firstvalue(k)= A(5);        % first integer value of signal (to test for errors)
end;
fclose(fid1);
clear A fid1;
%
%------ CARGA DATOS BINARIOS ----------------------------------------------
%SEÑAL 1
%-------
%-------ECG
%-------

signalECG= fullfile(PATH, DATAFILEECG);            % data in format 212
fid2=fopen(signalECG,'r');                         %'r' open file for reading
A= fread(fid2, [2,inf], 'uint8')';                 % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);

M1HA= bitand(A(:,2), 15);
PRLA=bitshift(bitand(A(:,2),8),9);     % sign-bit
ECG( : , 1)= bitshift(M1HA,8)+ A(:,1)-PRLA;
if ECG(1,1) ~= firstvalue, error('inconsistency in the first bit values'); end;

clear M1HA PRLA A;

%-------
%------- RES
%-------

signalRES= fullfile(PATH, DATAFILERES);            % data in format 16
fid3=fopen(signalRES,'r');                         %'r' open file for reading
B = fread(fid3,[4,inf], 'uint16')';                % matrix with 4 rows
fclose(fid3);
n=length(B);

RES= zeros(samples, 4);                                  % preallocating 

for k=1:n
    val1 = B(k,1);
    val2 = B(k,2);
    val3 = B(k,3);
    val4 = B(k,4);
    y1 = sign(2^(16-1)-val1)*(2^(16-1)-abs(2^(16-1)-val1)); 
    y2 = sign(2^(16-1)-val2)*(2^(16-1)-abs(2^(16-1)-val2)); 
    y3 = sign(2^(16-1)-val3)*(2^(16-1)-abs(2^(16-1)-val3)); 
    y4 = sign(2^(16-1)-val4)*(2^(16-1)-abs(2^(16-1)-val4)); 
    
    if ((y1 == 0) && (val1 ~= 0)) 
    RES(k , 1) = -val1; 
    else 
    RES(k , 1) = y1; 
    end
    
    if ((y2 == 0) && (val2 ~= 0)) 
    RES(k , 2) = -val2; 
    else 
    RES(k , 2) = y2; 
    end
    
    if ((y3 == 0) && (val3 ~= 0)) 
    RES(k , 3) = -val3; 
    else 
    RES(k , 3) = y3; 
    end
    
    if ((y4 == 0) && (val4 ~= 0))
    RES(k , 4) = -val4; 
    else 
    RES(k , 4) = y4; 
    end
end;

clear val1 val2 val3 val4 y1 y2 y3 y4 n fid2 fid3;

for k=2:nosig
    
    RES( : , k -1 )= (RES( : , k -1)- zerovalue(k))/gain(k);

end;
TIME=((0:(samples-1))/sfreq)';
signal1 =[TIME , ECG , RES];

clear ans bitres B dformat ECG firstvalue gain k nosig RES samples;
clear signalECG signalh signalRES z zerovalue;


fprintf(1,'\\n$> LOADING DATA FINISHED \n');

%------ CARGA LAS ANOTACIONES ---------------------------------------------
 
atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
fid3=fopen(atrd,'r');
A= fread(fid3, [2, inf], 'uint8')';
fclose(fid3);

ANNOT=[];
ATRTIME = [];
sa=size(A);
saa=sa(1);

i=1;
while i <= saa
    annoth = bitshift(A(i,2),-2);
    if annoth == 59
        ANNOT = [ANNOT;bitshift(A(i + 3,2),-2)];
      ATRTIME = [ATRTIME;A(i+2,1) + bitshift(A(i + 2,2),8) + bitshift(A(i + 1,1),16) + bitshift(A(i + 1,2),24)];
        i = i + 3;
    elseif annoth == 60
    elseif annoth == 61
    elseif annoth == 62
    elseif annoth == 63
        hilfe = bitshift(bitand(A(i,2),3),8) + A(i,1);
        hilfe = hilfe + mod(hilfe,2);
        i = i + hilfe/2;
    else
        ATRTIME = [ATRTIME;bitshift(bitand(A(i,2),3),8) + A(i,1)];
        ANNOT = [ANNOT;bitshift(A(i,2),-2)];
   end;
   i = i + 1;
end;

ANNOT(length(ANNOT)) = [];                  % Last Line = EOF (= 0)
ATRTIME(length(ATRTIME)) = [];              % Last Line = EOF
ATRTIME = (cumsum(ATRTIME))/sfreq;
ind = find(ATRTIME <= TIME(end));
ATRTIMED = ATRTIME(ind);
ANNOT = round(ANNOT);
ANNOTD = ANNOT(ind);
apnea1 = [ATRTIMED , ANNOTD];

clear DATAFILECG %, DATAFILERES, ARTFILE, HEADERFILE;
clear TIME;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%SEÑAL 2
% 
% sig='b01';
% 
% DATAFILEECG=[sig '.dat'];         % data-file
% DATAFILERES=[sig 'r' '.dat'];         % data-file
% ATRFILE = [sig 'r' '.apn'];       % apnea annotations
% HEADERFILE= [sig 'er' '.hea'];      % header-file in text format
% 
% %------ CARGA DATOS BINARIOS ----------------------------------------------
% % ------ CARGA ENCABEZADO DE LOS DATOS ------------------------------------
% 
% fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE);
% signalh= fullfile(PATH, HEADERFILE);
% fid1=fopen(signalh,'r');
% z= fgetl(fid1);
% A= sscanf(z, '%*s %d %d %d',[1,3]);
% nosig= A(1);  % number of signals
% sfreq=A(2);   % sample rate of data
% samples = A(3);
% clear A;
% for k=1:nosig
%     z= fgetl(fid1);
%     A= sscanf(z, '%*s %d %d %d %d %d %d %d %*s',[1,5]);
%     dformat(k)= A(1);           % format; 
%     gain(k)= A(2);              % number of integers per mV
%     bitres(k)= A(3);            % bitresolution
%     zerovalue(k)= A(4);         % integer value of ECG zero point
%     firstvalue(k)= A(5);        % first integer value of signal (to test for errors)
% end;
% fclose(fid1);
% clear A fid1;
% %
% %------ CARGA DATOS BINARIOS ----------------------------------------------
% %SEÑAL 2
% %-------
% %-------ECG
% %-------
% 
% signalECG= fullfile(PATH, DATAFILEECG);            % data in format 212
% fid2=fopen(signalECG,'r');                         %'r' open file for reading
% A= fread(fid2, [2,inf], 'uint8')';                 % matrix with 3 rows, each 8 bits long, = 2*12bit
% fclose(fid2);
% 
% M1HA= bitand(A(:,2), 15);
% PRLA=bitshift(bitand(A(:,2),8),9);     % sign-bit
% ECG( : , 1)= bitshift(M1HA,8)+ A(:,1)-PRLA;
% if ECG(1,1) ~= firstvalue, error('inconsistency in the first bit values'); end;
% 
% clear M1HA PRLA A;
% 
% %-------
% %------- RES
% %-------
% 
% signalRES= fullfile(PATH, DATAFILERES);            % data in format 16
% fid3=fopen(signalRES,'r');                         %'r' open file for reading
% B = fread(fid3,[4,inf], 'uint16')';                % matrix with 4 rows
% fclose(fid3);
% n=length(B);
% 
% RES= zeros(samples, 4);                                  % preallocating 
% 
% for k=1:n
%     val1 = B(k,1);
%     val2 = B(k,2);
%     val3 = B(k,3);
%     val4 = B(k,4);
%     y1 = sign(2^(16-1)-val1)*(2^(16-1)-abs(2^(16-1)-val1)); 
%     y2 = sign(2^(16-1)-val2)*(2^(16-1)-abs(2^(16-1)-val2)); 
%     y3 = sign(2^(16-1)-val3)*(2^(16-1)-abs(2^(16-1)-val3)); 
%     y4 = sign(2^(16-1)-val4)*(2^(16-1)-abs(2^(16-1)-val4)); 
%     
%     if ((y1 == 0) && (val1 ~= 0)) 
%     RES(k , 1) = -val1; 
%     else 
%     RES(k , 1) = y1; 
%     end
%     
%     if ((y2 == 0) && (val2 ~= 0)) 
%     RES(k , 2) = -val2; 
%     else 
%     RES(k , 2) = y2; 
%     end
%     
%     if ((y3 == 0) && (val3 ~= 0)) 
%     RES(k , 3) = -val3; 
%     else 
%     RES(k , 3) = y3; 
%     end
%     
%     if ((y4 == 0) && (val4 ~= 0))
%     RES(k , 4) = -val4; 
%     else 
%     RES(k , 4) = y4; 
%     end
% end;
% 
% clear val1 val2 val3 val4 y1 y2 y3 y4 n fid2 fid3;
% 
% for k=2:nosig
%     
%     RES( : , k -1 )= (RES( : , k -1)- zerovalue(k))/gain(k);
% 
% end;
% TIME=((0:(samples-1))/sfreq)';
% signal2 =[TIME , ECG , RES];
% 
% clear ans bitres B dformat ECG firstvalue gain k nosig RES samples;
% clear signalECG signalh signalRES z zerovalue;
% 
% 
% fprintf(1,'\\n$> LOADING DATA FINISHED \n');
% 
% %------ CARGA LAS ANOTACIONES ---------------------------------------------
%  
% atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
% fid3=fopen(atrd,'r');
% A= fread(fid3, [2, inf], 'uint8')';
% fclose(fid3);
% 
% ANNOT=[];
% ATRTIME = [];
% sa=size(A);
% saa=sa(1);
% 
% i=1;
% while i <= saa
%     annoth = bitshift(A(i,2),-2);
%     if annoth == 59
%         ANNOT = [ANNOT;bitshift(A(i + 3,2),-2)];
%       ATRTIME = [ATRTIME;A(i+2,1) + bitshift(A(i + 2,2),8) + bitshift(A(i + 1,1),16) + bitshift(A(i + 1,2),24)];
%         i = i + 3;
%     elseif annoth == 60
%     elseif annoth == 61
%     elseif annoth == 62
%     elseif annoth == 63
%         hilfe = bitshift(bitand(A(i,2),3),8) + A(i,1);
%         hilfe = hilfe + mod(hilfe,2);
%         i = i + hilfe/2;
%     else
%         ATRTIME = [ATRTIME;bitshift(bitand(A(i,2),3),8) + A(i,1)];
%         ANNOT = [ANNOT;bitshift(A(i,2),-2)];
%    end;
%    i = i + 1;
% end;
% 
% ANNOT(length(ANNOT)) = [];                  % Last Line = EOF (= 0)
% ATRTIME(length(ATRTIME)) = [];              % Last Line = EOF
% ATRTIME = (cumsum(ATRTIME))/sfreq;
% ind = find(ATRTIME <= TIME(end));
% ATRTIMED = ATRTIME(ind);
% ANNOT = round(ANNOT);
% ANNOTD = ANNOT(ind);
% apnea2 = [ATRTIMED , ANNOTD];
% 
% clear TIME 
% clear A ANNOT ANNOTD annoth ans atrd fid3 i ind sa;
% clear saa sig
% clear ATRFILE ATRTIME ATRTIMED DATAFILEECG DATAFILERES HEADERFILE
% 
% 
% 


%% ---------SEGMENTACIÓN -------------------
% ---------OBTENCIÓN DE CANDIDATOS
% ----- candidatos respecto a las anotaciones
% ---- señal 1

A1 = signal1(:,6);
A2 = signal1(:,2);

k=1;
for i=1:5999:length(A1)
    if i+5999<=length(A1)
        if k<= length(apnea1)
            spO2_1(:,k)=A1(i:i+5999)';
            ecg_1(:,k)=A2(i:i+5999)';
            k=k+1;
        end
    end
end

clear A1 A2;

% -----------señal 2
%  A1 = signal2(:,6);
%  A2 = signal2(:,2);
% k=1;
% for i=1:5999:length(A1)
%     if i+5999<=length(A1)
%         if k<= length(apnea2)
%             spO2_2(:,k)=A1(i:i+5999)';
%             ecg_2(:,k)=A2(i:i+5999)';
%             k=k+1;
%         end
%     end
% end


%unir las dos señales
%  apnea = [apnea1 ; apnea2];
%  ecg = [ecg_1 , ecg_2];  % a estas señales en particular a ecg_2 le falta un candidato
%  spO2 = [spO2_1 , spO2_2];

apnea = [apnea1];
ecg = [ecg_1];
spO2 = [spO2_1];

%clear A1 A2 i signal1 signal2 apnea1 apnea2 k spO2_1 spO2_2 ecg_1 ecg_2; 
% A1 = signal(:,6);
% A2 = signal(:,2);
% 
% k=1;
% for i=1:6000:length(A1)
%     if i+6000<=length(A1)
%         sp02(:,k)=A1(i:i+6000)';
%         ecg(:,k)=A2(i:i+6000)';
%         k=k+1;
%     end
% end

fprintf(1,'\\n$> SEGMENTATION DATA FINISHED \n');


%% --------------FEATURES ---------------------
% ...............SPO2 FEATURES
% ---------- The mean SpO2 value --------
M = mean (spO2);

% ---------- minimum SpO2 value ---------
MIN = min(spO2);

% % --------- the number of SpO2 values of less than
% % --------- 92 per ceent saturation ---------------

C = spO2 < 92;
LESS = sum(C);

%  % ................ECG FEATURES 
%  % --------------- PSD ecg -------------
% 
% %-------------------usando la función pwelch
 u=size(ecg);
 uu=u(2);
% 
n=length(apnea);

for i=1:uu
  [Pxx_u(:,i),f_u(:,i)]=pwelch(ecg(:,i),sfreq);
end
clear i;
% % % -------------DIFERENCIAS PSD 
% % %apnea
% % figure(1), clf, box on, hold on, grid on;
% % for i=100:800
% %   if apnea(i,2) == 8
% %     plot(f_u(:,1), pow2db(Pxx_u(:,i)));
% %     xlabel('Hz'); 
% %     ylabel('dB/Hz');
% %   end
% % end
% % %respiración normal
% % figure(2), clf, box on, hold on, grid on;
% % for i=1:uu
% %   if apnea(i,2) == 1
% %     plot(f_u(:,1), pow2db(Pxx_u(:,i)));
% %     xlabel('Hz'); 
% %     ylabel('dB/Hz');
% %   end
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
indVLF=find((f_u(:,1) >= 0.01) & (f_u(:,1) <= 0.045));
 
%hallar la energía en la banda de frecuencia VLF (0.019-0.036)
for i=1:uu
   enerVLF(i)=sum(Pxx_u(indVLF,i));
end
% 
% 
% %-------------------usando FFT
% % ecg = ecg'; %cambiar columnas por filas
% % ECG=fft(ecg);
% %  % como la señal fue muestreada a sfreq, entonces los indices
% %  % que representan los valores de Fourier deben escalarse para llevarlos
% %  % a escala de frecuencia
% %  
% %  N = length(ecg);
% %  frec = (0:N-1)*sfreq/N;
% % 
% %  % Hallar la PSD con el espectro completo (los dos lados simétricos)
% % Pxx =(abs(ECG).^2)/(N*sfreq);
% % if rem(N,2), % se evalua la paridad de N
% % select = 1:(N+1)/2; % si N es impar
% % else
% % select = 1:N/2+1; % si N es par
% % end
% % 
% %  %Para hallar la densidad de potencia desde 0 hasta la frecuencia de
% %  %Nyquist
% % Pxx=Pxx';
% % 
% % u=size(Pxx);
% % uu=u(2);
% % 
% % Pxx_unlado= zeros(length(select), uu);                                  % preallocating 
% % Pxx_u = zeros(length(select),uu);                                      %preallocating
% % for i=1:uu
% %     Pxx_unlado(:,i)=Pxx(select,i);  
% %     Pxx_u(:,i) = [Pxx_unlado(1,i); 2*Pxx_unlado(2:end-1,i); Pxx_unlado(end,i)];
% % end
% % clear i;
% % 
% % f_u = frec(select);
% %%hallar la energía en la banda de frecuencia VLF (0.013-0.0375)
% %indVLF=find((f_u >= 0.013) & (f_u <= 0.0375));
% 
% %EnerVLF = zeros(1,uu);                                     %preallocating
% 
% %for i=1:uu
%  %   EnerVLF(i)=sum(Pxx_u(indVLF,i));
% %end
% %EnerVLF=EnerVLF';
%  
% %------------------------------------------
% 
% 
% % figure (1), clf, box on, hold on, grid on;
% % plot(f_u(:,1), pow2db(Pxx_u(:,1)));
% % plot(f_u(:,100), pow2db(Pxx_u(:,100)));
% % plot(f_u(:,900), pow2db(Pxx_u(:,900)));
% % 
% % xlabel('Frequencia (Hz)')
% % ylabel('PSD (dB/Hz)')
% 
% 
% 
% % ------------------   detectar picos RR ----------------------------
% 
% %ecg = ecg'; %cambiar filas por columnas
% 
% % ------------------- graficar para un candidato
% 
% % [RR, t1]=rr_detect(ecg(:,1),sfreq);
% % TIME=((0:(length(ecg)-1))/sfreq)';
% 
% % figure(2), clf, box on, hold on, grid on;
% % plot(TIME,ecg(:,1));
% % title('ECG Y RR');
% % xlabel('Tiempo (s)');
% % ylabel('Voltios');
% % plot(t1, RR/100,'xk');
% 
% %clear RR t1 TIME;
% %--------------------desviación estandar de cada RR
% %--------------------mean(RR)
% % ------------------- RMSSD
u=size(ecg);
uu=u(2);

std_RR = zeros(1,uu);                                      %preallocating
mean_RR = zeros(1,uu);                                     %preallocating
RMSSD = zeros(1,uu);                                     %preallocating
for i=1:uu
   [RRa, t1a]=rr_detect(ecg(:,i),sfreq);

   for k=1:length(RRa)-1
    
       mean_RRa(k)=t1a(k+1)-t1a(k);
       mean_RRb(k)=mean_RRa(k)^2;
       if k == length(RRa)-1
           mean_RR(i)= mean(mean_RRa,2); % Promedio de toda la fila
           std_RR(i)=std(mean_RRa,0,2);
           RMSSD(i)=sqrt((1/(length(t1a)-1))*sum(mean_RRb,2));
       end 
   end
% % ------------------ energía banda VLF RR -------------   
%     %remuestreo 4Hz
%  %   tRRb = 0:0.25:t1a(end); 
%  %   RRb = spline(t1a,RRa,tRRb);
% 
%   %[Pxx_u,f_u]=pwelch(RRb,4); 
% % indVLF=find((f_u >= 0.003) & (f_u <= 0.04));
%  %  enerVLF(i)=sum(Pxx_u(indVLF));
%     
% 
% %figure (1), clf, box on, hold on, grid on;
% %plot(f_u, pow2db(Pxx_u));
% 
% 
 end
% clear mean_RRa mean_RRb;
% %
% % 
% % figure (1), clf, box on, hold on, grid on;
% % plot(f_u, pow2db(Pxx_u));
% % plot(f_u(:,100), pow2db(Pxx_u(:,100)));
% % plot(f_u(:,900), pow2db(Pxx_u(:,900)));
% % 
% % xlabel('Frequencia (Hz)')
% % ylabel('PSD (dB/Hz)')
% 
%--------- Graficar espacio de caracteristicas de dos 
%---------- dimensiones (M , MIN) 
% 
% n=u(2);
% %n=length(ecg);
% y1=apnea(1:n,2);
% mean_spo2=M(1:n);
% min_spo2=MIN(1:n);
% less_spo2=LESS(1:n);
% enerVLF=enerVLF(1:n);
% mean_RR=mean_RR(1:n);
% std_RR=std_RR(1:n);
% RMSSD = RMSSD(1:n);

% % % ---------- dimensiones (M , MIN) 
% X1=[mean_spo2;min_spo2]';
% figure(3), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea','Location','northwest');
%  xlabel('Valor medio de SpO2 (%)'); 
%  ylabel('Minimo valor de SpO2 (%)');
% 

% % 
% % % ---------- dimensiones (M , LESS) 
% % % 
%  X1=[mean_spo2;less_spo2]';
%  figure(4),  clf, box on, hold on, grid on;
%  xlabel('Valor medio de SpO2 (%)'); 
%  ylabel('Numero de valores menores al 92 en 1 minuto');
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea','Location','northwest');
% 
%  % plot(X1(1,y1==1),X1(2,y1==1),'r+',X1(1,y1==8),X1(2,y1==8),'bo')
% 
% % ---------- dimensiones (MIN , LESS) 
%  
%  X1=[min_spo2;less_spo2]';
%  figure(5),  clf, box on, hold on, grid on;
%  xlabel('Minimo valor de SpO2 (%)'); 
%  ylabel('Numero de valores menores al 92 en 1 minuto');
% %legend('normal','apnea','Location','northest');
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  % plot(X1(1,y1==1),X1(2,y1==1),'r+',X1(1,y1==8),X1(2,y1==8),'bo')
% 
% 
% %-------------------------(std(RR), EnerVLF)
% 
% X1=[std_RR;enerVLF]';
%  
% figure(6), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Desviación típica de los picos RR (ms)'); 
%  ylabel('FMB (ms^2)');
% 
% %-------------------------(std(RR), mean(RR))
% X1=[std_RR; mean_RR]';
% figure(7), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Desviación típica de los picos RR (ms)'); 
%  ylabel('Valor medio de los picos RR (ms)');
% 
%  %-------------------------(EnerVLF, mean(RR))
% X1=[enerVLF; mean_RR]';
% figure(8), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('FMB (ms^2)'); 
%  ylabel('Valor medio de los picos RR (ms)');
% 
% 
% %-------------------------(min_spo2,stdRR)
% X1=[min_spo2; std_RR]';
% figure(9), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Minimo valor de SpO2 (%)'); 
%  ylabel('Desviación típica de los picos RR (ms)');
% 
% %-------------------------(min_spo2,enerVLF)
% X1=[min_spo2; enerVLF]';
% figure(10), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Minimo valor de SpO2 (%)'); 
%  ylabel('FMB (ms^2)');
% 
% %-------------------------(min_spo2,meanRR)
% X1=[min_spo2; mean_RR]';
% figure(11), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Minimo valor de SpO2 (%)'); 
%  ylabel('Valor medio de los picos RR (ms)');
% 
% %-------------------------(mean_spo2,stdRR)
% X1=[mean_spo2; std_RR]';
% figure(12), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Valor medio de SpO2 (%)'); 
%  ylabel('Desviación típica de los picos RR (ms)');
% 
% %-------------------------(mean_spo2,enerVLF)
% X1=[mean_spo2; enerVLF]';
% figure(13), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Valor medio de SpO2 (%)'); 
%  ylabel('FMB (ms^2)');
% 
% %-------------------------(mean_spo2,mean_RR)
% X1=[mean_spo2; mean_RR]';
% figure(14), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Valor medio de SpO2 (%)'); 
%  ylabel('Valor medio de los picos RR (ms)');
%  
% %-------------------------(less_spo2,std_RR)
% X1=[less_spo2; std_RR]';
% figure(15), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Numero de valores menores al 92 en 1 minuto'); 
%  ylabel('Desviación típica de los picos RR (ms)');
%  
% %-------------------------(less_spo2,enerVLF)
% X1=[less_spo2; enerVLF]';
% figure(16), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Numero de valores menores al 92 en 1 minuto'); 
%  ylabel('FMB (ms^2)');
%  
% %-------------------------(less_spo2,mean_RR)
% X1=[less_spo2; mean_RR]';
% figure(17), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('Numero de valores menores al 92 en 1 minuto'); 
%  ylabel('Valor medio de los picos RR (ms)');
%  
% %--------------------------(RMSSD,std_RR)
% X1=[RMSSD; std_RR]';
% figure(18), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('RMSSD (ms)'); 
%  ylabel('desviación típica RR (ms)');
% 
% %--------------------------(RMSSD,enerVLF)
% X1=[RMSSD; enerVLF]';
% figure(19), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('RMSSD (ms)'); 
%  ylabel('FMB (ms^2)');
%  
% %--------------------------(RMSSD,mean_RR)
% X1=[RMSSD; mean_RR]';
% figure(20), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('RMSSD (ms)'); 
%  ylabel('promedio picos RR (ms)');
%  
% 
% %--------------------------(RMSSD,less_spo2)
% X1=[RMSSD; less_spo2]';
% figure(21), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('RMSSD (ms)'); 
%  ylabel('Numero de valores menores al 92 en 1 minuto');
% 
% 
% %--------------------------(RMSSD,mean_spo2)
% X1=[RMSSD; mean_spo2]';
% figure(22), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('RMSSD (ms)'); 
%  ylabel('promedio SpO_2 (%)');
% 
%  
% %--------------------------(RMSSD,min_spo2)
% X1=[RMSSD; min_spo2]';
% figure(23), clf, box on, hold on, grid on;
% gscatter(X1(:,1),X1(:,2),y1,'br','xo');
% legend('normal','apnea');
%  xlabel('RMSSD (ms)'); 
%  ylabel('minimo valor SpO2 (%)');
%  
%  
% fprintf(1,'\\n$> FEATURES EXTRACTION FINISHED \n');
% 
% % % -------------VERIFICACIÓN DE NORMALIDAD -----
% % h_minspo2=lillietest(min_spo2');
% % 
% % [f,x_values] = ecdf(min_spo2');
% % J = plot(x_values,f);
% % hold on;
% % K = plot(x_values,normcdf(x_values),'r--');
% % set(J,'LineWidth',2);
% % set(K,'LineWidth',2);
% % legend([J K],'min_spo2','Standard Normal CDF','Location','SE');
% % 
% % h_meanspo2=lillietest(mean_spo2');
% % h_lessspo2=lillietest(less_spo2');
% % h_stdRR=lillietest(std_RR');
% % h_enerVLF=lillietest(enerVLF');
% % h_meanRR=lillietest(mean_RR');
% % h_RMSSD=lillietest(RMSSD');
% % 
% % 
% % fprintf(1,'\\n$> KOLMOGOROV-SMIRNOV TEST FINISHED \n');
% % 
% % 
% %% -------------DISCRIMINANTE LINEAL ----------------------------
% % % ---------- dimensiones (M , MIN) 
% X1=[mean_spo2;min_spo2]';
% 
% lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% ldaResubErr1 = resubLoss(lda);      % Resubstitution error 
% ldaClass = resubPredict(lda);
% [ldaResubCM1,grpOrder1] = confusionmat(y1,ldaClass); % confusion matrix
% [confusionmatres1] = confusionmatStats(y1,ldaClass); % Evaluacion
% 
% figure (3),
% % Plot the curve K + [x,y]*L  = 0.
% f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% h2=ezplot(f,[-10 110 0 90]);
% h2.Color = 'r';
% h2.LineWidth = 2;
%  xlabel('Valor medio de SpO2 (%)'); 
%  ylabel('Minimo valor de SpO2 (%)');
% title('{\bf Clasificación con LDA}');
% 
% % 
% % % % ---------- dimensiones (M , LESS) 
% % X1=[mean_spo2;less_spo2]';
% % 
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr2 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM2,grpOrder2] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres2] = confusionmatStats(y1,ldaClass);
% % 
% % figure (4),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 110 0 7000]);
% % h2.Color = 'r';
% % h2.LineWidth = 2;
% %  xlabel('Valor medio de SpO2 (%)'); 
% %  ylabel('Numero de valores menores al 92 en 1 minuto');
% %  title('{\bf Clasificación con LDA}')
% % % % ---------- dimensiones (MIN , LESS) 
% %  X1=[min_spo2;less_spo2]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr3 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM3,grpOrder3] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres3] = confusionmatStats(y1,ldaClass);
% % 
% % figure (5),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 100 0 7000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % xlabel('Minimo valor de SpO2 (%)'); 
% %  ylabel('Numero de valores menores al 92 en 1 minuto');
% %  title('{\bf Clasificación con LDA}');
% %  
% %  %-------------------------(std(RR), EnerVLF)
% % 
% % X1=[std_RR;enerVLF]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr4 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM4,grpOrder4] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres4] = confusionmatStats(y1,ldaClass);
% % 
% % figure (6),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-20 2000 0 7000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %  xlabel('Desviación típica de los picos RR (ms)'); 
% %  ylabel('FMB (ms²)');
% %  title('{\bf Clasificación con LDA}');
% % 
% % %-------------------------(std(RR), mean(RR))
% % X1=[std_RR; mean_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr5 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM5,grpOrder5] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres5] = confusionmatStats(y1,ldaClass);
% % 
% % figure (7),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 2000 0 4000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %  xlabel('Desviación típica de los picos RR (ms)'); 
% %  ylabel('Valor medio de los picos RR (ms)');
% %  title('{\bf Clasificación con LDA}');
% %  
% %   %-------------------------(EnerVLF, mean(RR))
% % X1=[enerVLF; mean_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr6 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM6,grpOrder6] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres6] = confusionmatStats(y1,ldaClass);
% % 
% % figure (8),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-2000 7000 0 7000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %  xlabel('FMB (ms²)'); 
% %  ylabel('Valor medio de los picos RR (ms)');
% %  title('{\bf Clasificación con LDA}');
% % 
% % %-------------------------(min_spo2,stdRR)
% % X1=[min_spo2; std_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr7 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM7,grpOrder7] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres7] = confusionmatStats(y1,ldaClass);
% % 
% % figure (9),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 100 0 2000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % xlabel('Minimo valor de SpO2 (%)'); 
% %  ylabel('Desviación típica de los picos RR (ms)');
% %  title('{\bf Clasificación con LDA}');
% % 
% %  %-------------------------(min_spo2,enerVLF)
% % X1=[min_spo2; enerVLF]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr8 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM8,grpOrder8] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres8] = confusionmatStats(y1,ldaClass);
% % 
% % figure (10),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 100 0 70000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %  xlabel('Minimo valor de SpO2 (%)'); 
% %  ylabel('FMB (ms²)');
% %  title('{\bf Clasificación con LDA}');
% % 
% %  %-------------------------(min_spo2,meanRR)
% % X1=[min_spo2; mean_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr9 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM9,grpOrder9] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres9] = confusionmatStats(y1,ldaClass);
% % 
% % figure (11),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 100 0 3500]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %  xlabel('Minimo valor de SpO2 (%)'); 
% %  ylabel('Valor medio de los picos RR (ms)');
% %  title('{\bf Clasificación con LDA}');
% % 
% %  %-------------------------(mean_spo2,stdRR)
% % X1=[mean_spo2; std_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr10 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM10,grpOrder10] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres10] = confusionmatStats(y1,ldaClass);
% % 
% % figure (12),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 100 0 2000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % xlabel('Valor medio de SpO2 (%)'); 
% % ylabel('Desviación típica de los picos RR (ms)');
% % title('{\bf Clasificación con LDA}');
% % %-------------------------(mean_spo2,enerVLF)
% % X1=[mean_spo2; enerVLF]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr11 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM11,grpOrder11] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres11] = confusionmatStats(y1,ldaClass);
% % 
% % figure (13),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 100 0 70000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % xlabel('Valor medio de SpO2 (%)'); 
% %  ylabel('FMB (ms²)');
% % title('{\bf Clasificación con LDA}');
% % %-------------------------(mean_spo2,mean_RR)
% % X1=[mean_spo2; mean_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr12 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM12,grpOrder12] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres12] = confusionmatStats(y1,ldaClass);
% % 
% % figure (14),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 100 0 3500]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %  xlabel('Valor medio de SpO2 (%)'); 
% %  ylabel('Valor medio de los picos RR (ms)');
% % title('{\bf Clasificación con LDA}');
% % 
% % %-------------------------(less_spo2,std_RR)
% % X1=[less_spo2; std_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr13 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM13,grpOrder13] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres13] = confusionmatStats(y1,ldaClass);
% % 
% % figure (15),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 6000 0 2000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % xlabel('Numero de valores menores al 92 en 1 minuto'); 
% %  ylabel('Desviación típica de los picos RR (ms)');
% %  title('{\bf Clasificación con LDA}');
% %  
% % %-------------------------(less_spo2,enerVLF)
% % X1=[less_spo2; enerVLF]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr14 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM14,grpOrder14] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres14] = confusionmatStats(y1,ldaClass);
% % 
% % figure (16),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 7000 0 70000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % xlabel('Numero de valores menores al 92 en 1 minuto'); 
% %  ylabel('FMB (ms²)');
% %   title('{\bf Clasificación con LDA}');
% % 
% % %-------------------------(less_spo2,mean_RR)
% % X1=[less_spo2; mean_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr15 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM15,grpOrder15] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres15] = confusionmatStats(y1,ldaClass);
% % 
% % figure (17),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 2000 0 3500]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % xlabel('Numero de valores menores al 92 en 1 minuto'); 
% %  ylabel('Valor medio de los picos RR (ms)');
% %    title('{\bf Clasificación con LDA}');
% % 
% % %--------------------------(RMSSD,std_RR)
% % X1=[RMSSD; std_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr16 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM16,grpOrder16] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres16] = confusionmatStats(y1,ldaClass);
% % 
% % figure (18),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 1500 0 2500]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % title('{\bf Clasificación con LDA}');
% % xlabel('RMSSD (ms)'); 
% %  ylabel('Desviación típica RR (ms)');
% %    title('{\bf Clasificación con LDA}');
% % 
% % %--------------------------(RMSSD,enerVLF)
% % X1=[RMSSD; enerVLF]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr17 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM17,grpOrder17] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres17] = confusionmatStats(y1,ldaClass);
% % 
% % figure (19),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 1500 0 3500]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % title('{\bf Clasificación con LDA}');
% % xlabel('RMSSD (ms)'); 
% %  ylabel('FMB (ms²)');
% %    title('{\bf Clasificación con LDA}');
% % 
% % %--------------------------(RMSSD,mean_RR)
% % X1=[RMSSD; mean_RR]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr18 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM18,grpOrder18] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres18] = confusionmatStats(y1,ldaClass);
% % 
% % figure (20),
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[-10 1500 0 3500]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% % title('{\bf Clasificación con LDA}');
% % xlabel('RMSSD (ms)'); 
% %  ylabel('Promedio intervalo RR (ms)');
% %    title('{\bf Clasificación con LDA}');
% % 
% %    %--------------------------(RMSSD,less_spo2)
% % X1=[RMSSD; less_spo2]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr19 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM19,grpOrder19] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres19] = confusionmatStats(y1,ldaClass);
% % 
% % figure(21), 
% % % Plot the curve K + [x,y]*L  = 0.
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 6 0 6000]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %    title('{\bf Clasificación con LDA}');
% % legend('normal','apnea');
% %  xlabel('RMSSD (ms)'); 
% %  ylabel('Numero de valores menores al 92 en 1 minuto');
% % 
% % 
% % %--------------------------(RMSSD,mean_spo2)
% % X1=[RMSSD; mean_spo2]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr20 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM20,grpOrder20] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres20] = confusionmatStats(y1,ldaClass);
% % 
% % figure(22), 
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 6 0 100]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %    title('{\bf Clasificación con LDA}');
% % legend('normal','apnea');
% %  xlabel('RMSSD (ms)'); 
% %  ylabel('promedio SpO_2 (%)');
% % 
% %  
% % %--------------------------(RMSSD,min_spo2)
% % X1=[RMSSD; min_spo2]';
% % lda = fitcdiscr(X1,y1);             % Se crea clasificador lineal
% % K= lda.Coeffs(1,2).Const;           % Entrega los coeficientes
% % L= lda.Coeffs(1,2).Linear;          % Limite entre las clases
% % ldaResubErr21 = resubLoss(lda);      % Resubstitution error 
% % ldaClass = resubPredict(lda);
% % [ldaResubCM21,grpOrder21] = confusionmat(y1,ldaClass); % confusion matrix
% % [confusionmatres21] = confusionmatStats(y1,ldaClass);
% % 
% % 
% % figure(23), 
% % f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
% % h2=ezplot(f,[0 6 0 100]);
% % h2.Color = 'k';
% % h2.LineWidth = 2;
% %    title('{\bf Clasificación con LDA}');
% % legend('normal','apnea');
% %  xlabel('RMSSD (ms)'); 
% %  ylabel('minimo valor SpO2 (%)');
% %  
% % %  
