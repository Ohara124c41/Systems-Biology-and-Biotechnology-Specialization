
colors = repmat('krgbmc',1,300) ;

Rtot = 20 ;

L = 0:0.01:200 ;

KDs = 10:20:90 ;

figure
hold on

for i=1:length(KDs) 
  KD = KDs(i) ;
  LR = Rtot*L./(L + KD) ;
  plot(L,LR,colors(i)) ;
  LR_all(i,:) = LR ;
  figurelegend{i} = ['K_D = ',int2str(KD),' uM'] ;
end
  
xlabel('[Ligand] (uM)')
ylabel('[Ligand-Receptor] (nM)')
legend(figurelegend,'Location','SouthEast')

