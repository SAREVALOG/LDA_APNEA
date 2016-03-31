function [RR t] = rr_pan_tompkin(ecg, fm);
%% function [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg,fs)
% Complete implementation of Pan-Tompkins algorithm

%% Inputs
% ecg : raw ecg vector signal 1d signal
% fs : sampling frequency e.g. 200Hz, 400Hz and etc
% gr : flag to plot or not plot (set it 1 to have a plot or set it zero not
% to see any plots
%% Outputs
% qrs_amp_raw : amplitude of R waves amplitudes
% qrs_i_raw : index of R waves
% delay : number of samples which the signal is delayed due to the
% filtering
%% Method :

%% PreProcessing
% 1) Signal is preprocessed , if the sampling frequency is higher then it is downsampled
% and if it is lower upsampled to make the sampling frequency 200 Hz
% with the same filtering setups introduced in Pan
% tompkins paper (a combination of low pass and high pass filter 5-15 Hz)
% to get rid of the baseline wander and muscle noise. 

% 2) The filtered signal
% is derivated using a derivating filter to high light the QRS complex.

% 3) Signal is squared.4)Signal is averaged with a moving window to get rid
% of noise (0.150 seconds length).

% 5) depending on the sampling frequency of your signal the filtering
% options are changed to best match the characteristics of your ecg signal

% 6) Unlike the other implementations in this implementation the desicion
% rule of the Pan tompkins is implemented completely.

%% Decision Rule 
% At this point in the algorithm, the preceding stages have produced a roughly pulse-shaped
% waveform at the output of the MWI . The determination as to whether this pulse
% corresponds to a QRS complex (as opposed to a high-sloped T-wave or a noise artefact) is
% performed with an adaptive thresholding operation and other decision
% rules outlined below;

% a) FIDUCIAL MARK - The waveform is first processed to produce a set of weighted unit
% samples at the location of the MWI maxima. This is done in order to localize the QRS
% complex to a single instant of time. The w[k] weighting is the maxima value.

% b) THRESHOLDING - When analyzing the amplitude of the MWI output, the algorithm uses
% two threshold values (THR_SIG and THR_NOISE, appropriately initialized during a brief
% 2 second training phase) that continuously adapt to changing ECG signal quality. The
% first pass through y[n] uses these thresholds to classify the each non-zero sample
% (CURRENTPEAK) as either signal or noise:
% If CURRENTPEAK > THR_SIG, that location is identified as a QRS complex
% candidate and the signal level (SIG_LEV) is updated:
% SIG _ LEV = 0.125 CURRENTPEAK + 0.875 SIG _ LEV

% If THR_NOISE < CURRENTPEAK < THR_SIG, then that location is identified as a
% noise peak and the noise level (NOISE_LEV) is updated:
% NOISE _ LEV = 0.125CURRENTPEAK + 0.875 NOISE _ LEV
% Based on new estimates of the signal and noise levels (SIG_LEV and NOISE_LEV,
% respectively) at that point in the ECG, the thresholds are adjusted as follows:
% THR _ SIG = NOISE _ LEV + 0.25  (SIG _ LEV ? NOISE _ LEV )
% THR _ NOISE = 0.5 (THR _ SIG)
% These adjustments lower the threshold gradually in signal segments that are deemed to
% be of poorer quality.


% c) SEARCHBACK FOR MISSED QRS COMPLEXES - In the thresholding step above, if
% CURRENTPEAK < THR_SIG, the peak is deemed not to have resulted from a QRS
% complex. If however, an unreasonably long period has expired without an abovethreshold
% peak, the algorithm will assume a QRS has been missed and perform a
% searchback. This limits the number of false negatives. The minimum time used to trigger
% a searchback is 1.66 times the current R peak to R peak time period (called the RR
% interval). This value has a physiological origin - the time value between adjacent
% heartbeats cannot change more quickly than this. The missed QRS complex is assumed
% to occur at the location of the highest peak in the interval that lies between THR_SIG and
% THR_NOISE. In this algorithm, two average RR intervals are stored,the first RR interval is 
% calculated as an average of the last eight QRS locations in order to adapt to changing heart 
% rate and the second RR interval mean is the mean 
% of the most regular RR intervals . The threshold is lowered if the heart rate is not regular 
% to improve detection.

% d) ELIMINATION OF MULTIPLE DETECTIONS WITHIN REFRACTORY PERIOD - It is
% impossible for a legitimate QRS complex to occur if it lies within 200ms after a previously
% detected one. This constraint is a physiological one  due to the refractory period during
% which ventricular depolarization cannot occur despite a stimulus[1]. As QRS complex
% candidates are generated, the algorithm eliminates such physically impossible events,
% thereby reducing false positives.

% e) T WAVE DISCRIMINATION - Finally, if a QRS candidate occurs after the 200ms
% refractory period but within 360ms of the previous QRS, the algorithm determines
% whether this is a genuine QRS complex of the next heartbeat or an abnormally prominent
% T wave. This decision is based on the mean slope of the waveform at that position. A slope of
% less than one half that of the previous QRS complex is consistent with the slower
% changing behaviour of a T wave  otherwise, it becomes a QRS detection.
% Extra concept : beside the points mentioned in the paper, this code also
% checks if the occured peak which is less than 360 msec latency has also a
% latency less than 0,5*mean_RR if yes this is counted as noise

% f) In the final stage , the output of R waves detected in smoothed signal is analyzed and double
% checked with the help of the output of the bandpass signal to improve
% detection and find the original index of the real R waves on the raw ecg
% signal

%% References :

%[1]PAN.J, TOMPKINS. W.J,"A Real-Time QRS Detection Algorithm" IEEE
%TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. BME-32, NO. 3, MARCH 1985.

%% Author : Hooman Sedghamiz
% Linkoping university 
% email : hoose792@student.liu.se
% hooman.sedghamiz@medel.com

% Any direct or indirect use of this code should be referenced 
% Copyright march 2014
%%
[idx_RR] = decgmf_curso(60*ecg);

RR = diff(abs(idx_RR)/fm)*1000;

t = cumsum(RR/1000);



function [r] = decgmf_curso(ecg,or_qrs,mu_qrs);
% [r,rm] = decgmf1(ecg,or_qrs,mu_qrs);
% Deteccion de onda R usando un filtro adaptable para delimitar una region de
% busqueda de la onda. La deteccion dentro de la region se realiza por maximo
% de la se~al original o por maximo de la se~al filtrada por templete 
% (promedio de ultimas 5 ondas QRS detectadas).
% r : deteccion por maximo
% rm : deteccion por filtrado por templete 
% ecg : se~al original
% or_qrs : orden del filtro adaptable
% mu_qrs : constante de convergencia

% parametros para la localizacion de la region QRS
if nargin < 3,			% si el numero de argumentos es menor que 3
 or_qrs = 32;
			%   orden de default
 mu_qrs = 0.0001;		%   constante de convergencia de default
end;				% fin de valores de default

% valores iniciales
fs = 200;    %%3*1024;			% todo se corrio con 200!
w_qrs = zeros(or_qrs+1,1);
	% pesos (+1 por termino de sesgo)
reg = zeros(size(ecg));		% qrs region 
	
auxu = zeros(or_qrs+1,1);
	% variable auxiliar: entrada al banco de filtros
auxe = zeros(size(ecg));
	% error del filtrado LMS
% parametros para la deteccion 
level0 = 200;		% umbral minimo para el inicio de la region de busqueda
level1 = 400;		% umbral inicial para la deteccion de la onda R
minPeriod = fs*0.25;	% minimo periodo (FM*60/240)
bandreg = 0;		% bandera: en la region de busqueda
bandr = 0;		% bandera: deteccion de onda R
RRok = 0;		% bandera: intervalo RR en rango 
band5 = 0;		% bandera: 5 latidos correctos
kreg = 0;		% indice al inicio de la region de deteccion
rr = 1;			% indice auxiliar para R
nr = 0;			% numero de ondas R detectada
load tw.dat		% templete original
wait = 0;		% bandera: espera
k = 0;			% contador
tf = zeros(size(ecg));	% salida del filtro por templete
template = [tw tw tw tw tw];	% arreglo para mantener ultimos cinco complejos
qrsi = 1;		% variable auxiliar de inicio de region
it = 1;			% indice al templete a remplazar en "template"
RR = ones(5,1);         % arreglo para mantener ultimos cinco intervalos RR
rrmean = 1;		% RR promedio de los ultimos cinco RRs
rrstd = 0;		% desviacion estandard de los ultimos cinco RRs
ampR = max(tw);		% maxima amplitud de region de busqueda de complejo
r(1) = 1;
rm(1) = 1;

% derivando el ECG
 decg = ecg(2:length(ecg))-ecg(1:length(ecg)-1);
 decg(length(ecg)) = decg(length(decg));
	% completa arreglo de derivada
 for i = or_qrs:length(ecg)
	% para cada dato

% Delimitacion de region por LMS  
 j = fix(i-or_qrs/2);
		%   indice interno
 auxu = [decg(i-or_qrs+1:i);1];	%   copia entrada (1 por termino de sesgo)
 reg(j) = auxu'*w_qrs; 		%   salida del filtro
 auxe = decg(j)-reg(j);		%   error
 w_qrs = w_qrs+2*mu_qrs*auxe*auxu;	% actualiza pesos del filtro
 s = sqrt(w_qrs'*w_qrs);	%   normalizacion de los pesos
 if s ~= 0 
  w_qrs = w_qrs/s;
 end; 
 if i > or_qrs
			%   define region de busqueda  
  j = j-1;
  reg(j) = abs(reg(j+1) - reg(j)).^2/100;	% como c*abs(D(salida filtro))^2
% deteccion del complejo QRS 
  lt = length(tw);		%   longitud del templete 
  trend = sum(reg(j-4:j) >= level0)/5;	% fraccion de puntos reg(j-4:j)>level0
  if ((trend == 1&bandreg == 0&j-5-qrsi > minPeriod)|...
      (reg(j)>=level1&bandreg == 1)|(wait == 1))&(j-5+lt/2 < length(ecg))     %   si ((5 puntos consecutivos > level0 y 
      				%    bandera de region=0  y periodo minimo 
      				%    se cumple) o (reg > umbral de deteccion
      				%    de deteccion de R o espera = 1)) y no en el
      				%    ultimo segmento del archivo 	   
   kreg = j-5;	 		%     indice al inicio de la region  
   if bandreg == 0		%     si entrando a la region
    level1 = 0.2*reg(j);
    nr = nr + 1;
		%       incrementa contador de latidos
    qrsi = kreg;
		%       copia inicio a variable auxiliar
    bandreg = 1;
		%       enciende bandera de region
   end;
 			%     fin si entrando a region
   wait = 1;
			%     espera = 1
%    if(kreg<lt/2)
% 		%     filtrado por templete
%     tf(kreg) = tw'*[zeros(fix(lt/2-kreg),1);ecg(fix(1:kreg+lt/2))];
%  	
%    else 
%  
  %  tf(kreg) = tw'*(ecg(fix(kreg-lt/2+1:kreg+lt/2))); 
  %end;
			%     fin de filtrado
   if (reg(j) < level1 & bandreg == 1 & trend < 1)
   
   				%     si reg < umbral de deteccion de R y ya se
   				%      entro a la region y 5 puntos < level0
   				%      (posible salida de la region de busqueda)
    if k < lt/4
		%       si falsa salida en primer cuarto
     k = k + 1;
		%         incrementa contador
    else
			%       otro		
     bandreg = 0;
		%         salio de region de busqueda
     k = 0;
			%         inicializa contador
     wait = 0;
			%         espera = 0
     mqrs = mean(ecg(qrsi:kreg));	% calcula media del ECG en region
     [mm,pm] = max(reg(qrsi:kreg));
     if band5
      RRok =  (abs(abs(r(nr))-abs(r(nr-1))-rrmean) < rrstd);
     end;
     if mm > 0.4*ampR 		% era 0.6
      if mm < 1.4*ampR
       ampR = (ampR+mm)/2;
      end;
     if (band5 & RRok)|...
      (~band5 & max(ecg(qrsi:kreg)-mqrs) > abs(min(ecg(qrsi:kreg))-mqrs))
                                %         si (ya hubo mas de 5 latidos y nuevo
                                %          RR en media +/- sd) o (menos de 5
                                %          latidos y amplitud R > abs(QRSmin)				
        [m,r(nr)] = max(ecg(qrsi:kreg));	%   detecta maximo en region (R por max)
  %      [mf,rm(nr)] = max(tf(qrsi:kreg));	%   detecta maximo en filtro (R por ft)
        r(nr) = ceil(r(nr) + qrsi - 1);	%   calcula posicion de onda R (max)
  %      rm(nr) = ceil(rm(nr) + qrsi - 1); %   (ft = filtrado por templete)
        if nr == 2			%   actualiza estadisticos RR
         RR = (abs(r(nr)) - abs(r(nr-1)))*ones(5,1);
        elseif nr > 1
         RR(it) = abs(r(nr)) - abs(r(nr - 1));	 
        end; 
        rrmean = mean(RR);
        rrstd = std(RR);
        template(:,it) = ecg(fix(r(nr)-lt/2+1:r(nr)+lt/2));	% actualiza templete
        tw = mean(template')';
        it = rem(it,5) + 1;
        if band5 == 0 & it == 5
         band5 == 1;
        end;		
     else				%   posible foco ectopico
      [m,r(nr)] = max(ecg(qrsi:kreg)-mean(ecg(qrsi:kreg))); % detecta posible Rs
%      [m,rm(nr)] = max(abs(tf(qrsi:kreg)));	% pero no incluye valores para
      r(nr) = -ceil(r(nr) + qrsi - 1);		% templete o estadisticos RR
 %     rm(nr) = -ceil(rm(nr) + qrsi - 1);	% almacena las Rs con signo -
     end;			%         fin complejo ok
    else 
     nr = nr-1;
     qrsi = r(nr);
    end;
    end;
 			%       fin falsa salida 
   end;
			%     fin posible salida
  end;
				%   fin region de busqueda

 end;
				% fin i > orqrs
end; 
				% fin iteraciones
                
clear reg;
r = r(:);			% vectores columna
%rm = rm(:);
