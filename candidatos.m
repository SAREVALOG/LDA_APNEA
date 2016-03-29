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
% Extrae caracteristicas de los registros
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

%%%%%%%%%SEÑAL 2

sig='b01';

DATAFILEECG=[sig '.dat'];         % data-file
DATAFILERES=[sig 'r' '.dat'];         % data-file
ATRFILE = [sig 'r' '.apn'];       % apnea annotations
HEADERFILE= [sig 'er' '.hea'];      % header-file in text format

%------ CARGA DATOS BINARIOS ----------------------------------------------
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
%SEÑAL 2
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
signal2 =[TIME , ECG , RES];

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
apnea2 = [ATRTIMED , ANNOTD];

clear TIME 
clear A ANNOT ANNOTD annoth ans atrd fid3 i ind sa;
clear saa sig
clear ATRFILE ATRTIME ATRTIMED DATAFILEECG DATAFILERES HEADERFILE





% ---------SEGMENTACIÓN -------------------
%---------OBTENCIÓN DE CANDIDATOS
%----- candidatos respecto a las anotaciones
%---- señal 1

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

%-----------señal 2
 A1 = signal2(:,6);
 A2 = signal2(:,2);
k=1;
for i=1:5999:length(A1)
    if i+5999<=length(A1)
        if k<= length(apnea2)
            spO2_2(:,k)=A1(i:i+5999)';
            ecg_2(:,k)=A2(i:i+5999)';
            k=k+1;
        end
    end
end


%------------unir las dos señales
 apnea = [apnea1 ; apnea2];
 ecg = [ecg_1 , ecg_2];  % a estas señales en particular a ecg_2 le falta un candidato
 spO2 = [spO2_1 , spO2_2];

apnea = [apnea1];
ecg = [ecg_1];
spO2 = [spO2_1];

clear A1 A2 i signal1 signal2 apnea1 apnea2 k spO2_1 spO2_2 ecg_1 ecg_2; 
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
   [RRa, t1a]=rr_pan_tompkin(ecg(:,i),sfreq);

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
clear mean_RRa mean_RRb;
%
% 
% % figure (1), clf, box on, hold on, grid on;
% % plot(f_u, pow2db(Pxx_u));
% % plot(f_u(:,100), pow2db(Pxx_u(:,100)));
% % plot(f_u(:,900), pow2db(Pxx_u(:,900)));
% % 
% % xlabel('Frequencia (Hz)')
% % ylabel('PSD (dB/Hz)')
% 
