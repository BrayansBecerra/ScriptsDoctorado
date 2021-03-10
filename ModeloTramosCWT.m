%% definicion de modelo estacionario por tramos %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs=4;
 f1=0.003;
 f2=0.1;
 f3=0.2;
 
 t = 0:1/fs:10800;
 x = 0.2*sin(2*pi*f1*t)+0.15*sin(2*pi*f2*t)+ 0.1*sin(2*pi*f3*t); %%SEÑAL CONJUNTA DE 3 COMPONENTES

 figure(1)
 plot(t,x,'red')
%%%%%%%%%%%%%%%%%%%%%%% definicion de señal por tramos %%%%%%%%%%%%%%%%%%%

% %  fs=4;
% %  f1=0.003;
% %  f2=0.1;
% %  f3=0.2;
% %  %t = 0:1/fs:10800;
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
 [S,F,P]=cwt(x,'amor',fs);
 s=abs(S);

ii=figure(2)
 surf(t/60,F,s,'FaceColor','interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');
 axis tight;
view(0,90);
ylim([0 0.4]);
xlim([0 180]);

zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (min)');
title('Representación T-F de modelo estacionario con CWT');
colormap('jet');
colorbar
print(ii,'ScalModel3compCWT','-dtiff','-r300');

 