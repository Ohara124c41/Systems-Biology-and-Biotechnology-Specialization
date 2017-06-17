
colors = repmat('krgbmc',1,300) ;
% this is to define colors for plotting

% another way to do a comment

L = 0:0.01:200 ;
Rtot = 20 ;

figure
hold on 

KD = 10 ;
LR = Rtot*L./(L + KD) ;
plot(L,LR,colors(1))

KD = 30 ;
LR = Rtot*L./(L + KD) ;
plot(L,LR,colors(2))

KD = 50 ;
LR = Rtot*L./(L + KD) ;
plot(L,LR,colors(3))

KD = 70 ;
LR = Rtot*L./(L + KD) ;
plot(L,LR,colors(4))

KD = 90 ;
LR = Rtot*L./(L + KD) ;
plot(L,LR,colors(5))

xlabel('[Ligand] (uM)')
ylabel('[Ligand-Receptor] (nM)')
legend('K_D = 10 uM','K_D = 30 uM', ...
  'K_D = 50 uM','K_D = 70 uM', ...
  'K_D = 90 uM', ...
  'Location','SouthEast')

    



