# ***************************************************************************
#  coeficientes de friccion propios de la teoria F en combinacion con una 
#  ecuacion de estado cubica para el modelado de viscosidades dinamicas en
#  la fase liquida de liquidos ionicos (LI) y de solventes eutecticos 
#  profundos (DES).
#    Soave SRK (ede=1)
#    Peng-Robinson PR (ede=2)
#
#  Funciones Alfa consideradas:
#    Soave (alfa=1)
#    Stryjek-Vera (alfa=2)
#    Melhem et al. (alfa=3)
#    Yokozeki (alfa=4)
#
#  Parametros ajustables:
#    A1, A2, A3, B1, B2, B3, C2, C3 (poflag=1)
#    A1, A2, A3, B1, B2, B3, C4 (poflag=2)
#    A0, A1, A2, B0, B1, B2, C0, C1, C2 (poflag=3)
#    A1, A2, B1, B2, C2, C3 (poflag=4)
#    B1, C2, C3 (poflag=5)
#
#  Metodos de optimizacion considerados:
#    Levenberg-Marquardt (moflag=1)
#    Nelder-Mead (moflag=2)
#    Powell (moflag=3)
#    Differential-Evolution (moflag=4)
#
#  Funciones externas requeridas por EtaLIDES_TF:
#    prhoede, visco
#
c3 =  cell2mat(temp(4));
temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);
c4 =  cell2mat(temp(4));
temp = textscan(fid, '%s %d', 1,'HeaderLines', 7);
  etacal = visco(texp, x);
  y = etacal ./ etaexp;
      [fval, xfit, cvg, iter] = leasqr(texp, varobs, x, @flm);
          ctl.XVmax = [ 1  .1  .01  1  .1  .01  .01];
        case 3
          ctl.XVmax = [ 1  .1  .01  1  .1  .01  .1  .01  .001];
        case 4
    fprintf(fid, 'A2 [cP/bar]       %14.7e\n', xfit(2));
    fprintf(fid, 'A3 [cP/bar]       %14.7e\n', xfit(3));
    fprintf(fid, 'B1 [cP/bar]       %14.7e\n', xfit(4));
    fprintf(fid, 'B2 [cP/bar]       %14.7e\n', xfit(5));
    fprintf(fid, 'B3 [cP/bar]       %14.7e\n', xfit(6));
    fprintf(fid, 'C4 [cP/bar2]      %14.7e\n\n\n', xfit(7));
  case 3
    fprintf(fid, 'A1 [cP/bar]       %14.7e\n', xfit(2));
    fprintf(fid, 'A2 [cP/bar]       %14.7e\n', xfit(3));
    fprintf(fid, 'B0 [cP/bar]       %14.7e\n', xfit(4));
    fprintf(fid, 'B1 [cP/bar]       %14.7e\n', xfit(5));
    fprintf(fid, 'B2 [cP/bar]       %14.7e\n', xfit(6));
    fprintf(fid, 'C0 [cP/bar2]      %14.7e\n', xfit(7));
    fprintf(fid, 'C1 [cP/bar2]      %14.7e\n', xfit(8));
    fprintf(fid, 'C2 [cP/bar2]      %14.7e\n\n\n', xfit(9));
  case 4
  figure
  grid
  scatter3(pexp, texp, etaexp, 80, 'r', 'filled')
  xlabel('Pressure [bar]')
  ylabel('Temperature [K]')
  zlabel('Viscosity [mPa-s]')
  title(sprintf('Experimental & Calculated Viscosities vs T & P for %s', compuesto))
  hold on
  pp = linspace (min(pexp), max(pexp), 100)';
  surf(xx, yy, zz, 'FaceAlpha', .75, 'EdgeColor', 'none');
  set(gca,'ZScale','log')
  set(gcf, 'Position', [200, 200, 800, 800]);
  view(135, 25)
  set(gca, 'fontsize', 18);
  label1 = [model2, model1, ' Model'];
  legend('Experimental Data', label1)
  
  figure
  scatter(pexp, rdev_eta, 90, texp, 'filled')
  grid
  xlabel('Pressure [bar]')
  ylabel('Viscosity Rel. Dev. [%]')
  title(sprintf('Viscosity Relative Deviations vs P for %s, Colorbar represents T', compuesto))
  set(gcf, 'Position', [100, 100, 900, 600]);
  set(gca, 'fontsize', 18);
  c=colorbar();
  tmin = min(texp); tmax = max(texp);
  set(c, 'ytick', tmin:10:tmax)
  labels = {};
  for v = get(c, 'ytick')
    labels{end+1} = sprintf('%.0f K', v);
  endfor
  set(c, 'yticklabel', labels);
  set(c, 'fontsize', 18);
endif