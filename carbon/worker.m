% worker to call CO2SYS.m

dir0 = '../../ptools_output/carbon/temp/';
a = load([dir0,'input.mat']);


A = CO2SYS(a.alkalinity(:), a.TIC(:), 1, 2, a.salt(:), a.temp(:), a.temp(:), 0, 0, 50, 2, 1, 10, 1);

ph = A(:,18);
om = A(:,31);

PH = reshape(ph, size(a.salt));
OM = reshape(om, size(a.salt));

save([dir0,'PH.mat'], 'PH');
save([dir0,'OM.mat'], 'OM');