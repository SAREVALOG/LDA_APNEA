function [RR t] = rr_pan_tompkin(ecg, fm);

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
