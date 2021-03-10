%% definicion de modelo estacionario por tramos %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 fs=4;
 f1=0.003;
 f2=0.1;
 f3=0.2;
 
 t = 0:1/fs:10800;
 x = 0.2*sin(2*pi*f1*t)+0.15*sin(2*pi*f2*t)+ 0.1*sin(2*pi*f3*t); %%SEÑAL CONJUNTA DE 3 COMPONENTES
 
 %%%%%%%%%%%%%%%%%%%%%% SEGMENTO PARA EVALUACION DE SEÑAL POR TRAMOS %%%%%%
 
% %  t1=0:1/fs:(3600-(1/fs));
% %  t2=3600:1/fs:(7200-(1/fs));
% %  t3=7200:1/fs:(10800-(1/fs));
% %  
% %  x1 = 0.2*sin(2*pi*f1*t1);
% %  x2 = 0.15*sin(2*pi*f2*t2);
% %  x3 = 0.1*sin(2*pi*f3*t3);
% %  
% %  t=cat(1,t1',t2',t3');
% %  x=cat(1,x1',x2',x3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  


 kk=figure(1)
 plot(t/60,x,'red')
 xlabel('tiempo (min)');
 ylabel('Amplitud');
 xlim([0 180]);
 title('Señal de un modelo estacionario de 3 componentes')
 print(kk,'modelo3compo','-dtiff','-r300');
 
 
 
 
duracion=600; %% TAMAÑO DE LA VENTANA DE ANALISIS EN SEGUNDOS
ventana=duracion*fs;

  a=x(1:ventana/2);
  aa=flip(a);
  b=x(end-(ventana/2):end-1);
  c=flip(b);
  
X=cat(1,aa',x',c');
  
gridSegundos=360; %% DESPLAZAMIENTO DEL GRID EN SEGUNDOS
puntoscorri=(gridSegundos*4); %+((gridSegundos/2)*4)
grid=ventana-puntoscorri;
nfft=2*ventana; %% NUMERO DE PUNTOS DE LA STFT

 [S,F,T,P]=spectrogram(X,hann(ventana),grid,nfft,fs,'yaxis');
% % % 
ii=figure(2)
surf((T-T(1))/60,F,P,'FaceColor', 'interp',...
     'EdgeColor','none',...
     'FaceLighting','phong');
 axis tight;
 view(0,90);
 ylabel('Frecuencia (Hz)');
 xlabel('Tiempo (min)');
 title('Representación T-F de modelo estacionario con STFT');
 colorbar
 colormap('jet');
 ylim([0 0.4]);
% caxis([0 max(max(P(:,:)))/20]);

 print(ii,'modelo3compoSTFT','-dtiff','-r300');