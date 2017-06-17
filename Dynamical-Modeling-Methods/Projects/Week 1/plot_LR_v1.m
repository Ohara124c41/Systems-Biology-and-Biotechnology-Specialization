
L = 0:0.01:200 ;
Rtot = 20 ;

KD = 10 ;
LR = Rtot*L./(L + KD) ;
figure
plot(L,LR)
title('K_D = 10 uM')

KD = 30 ;
LR = Rtot*L./(L + KD) ;
figure
plot(L,LR)
title('K_D = 30 uM')

KD = 50 ;
LR = Rtot*L./(L + KD) ;
figure
plot(L,LR)
title('K_D = 50 uM')

KD = 70 ;
LR = Rtot*L./(L + KD) ;
figure
plot(L,LR)
title('K_D = 70 uM')

KD = 90 ;
LR = Rtot*L./(L + KD) ;
figure
plot(L,LR)
title('K_D = 90 uM')

    



