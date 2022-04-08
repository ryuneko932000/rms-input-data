# ***************************************************************************##                   PROGRAMA   E t a L I D E S _ T F##  EtaLIDES_TF es el programa principal para la optimizacion de los
#  coeficientes de friccion propios de la teoria F en combinacion con una 
#  ecuacion de estado cubica para el modelado de viscosidades dinamicas en
#  la fase liquida de liquidos ionicos (LI) y de solventes eutecticos 
#  profundos (DES).##  EdE cubicas consideradas:
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
## ***************************************************************************clear allclose allclcglobal pm tc pc omega kappa1 alfa0 alfa1 beta1 beta2 ctvglobal pexp texp etaexpglobal ede alfa poflagglobal pr padisp('**** Modelado de Viscosidade de LI & DES con la Teoria de Friccion ****') disp(' ')nombre_arch = input('Nombre del archivo de entrada => ', 's'); # leer el archivo de entradafid = fopen(nombre_arch);prop = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f', 1, 'HeaderLines', 1);compuesto = char(prop(1));pm = cell2mat(prop(2));tc = cell2mat(prop(3));pc = cell2mat(prop(4));omega = cell2mat(prop(5));kappa1 = cell2mat(prop(6));alfa0 = cell2mat(prop(7));alfa1 = cell2mat(prop(8));beta1 = cell2mat(prop(9));beta2 = cell2mat(prop(10));ctv =  cell2mat(prop(11));temp = textscan(fid, '%s %d', 1,'HeaderLines', 10);mtflag = cell2mat(temp(2));switch mtflag  case 1    ede = 1; alfa = 1; model1 = 'SRK-Soave';  case 2    ede = 2; alfa = 1; model1 = 'PR-Soave';  case 3    ede = 2; alfa = 2; model1 = 'PR-SV';  case 4    ede = 1; alfa = 3; model1 = 'SRK-Melhem';  case 5    ede = 2; alfa = 3; model1 = 'PR-Melhem';  case 6    ede = 1; alfa = 4; model1 = 'SRK-Yokozeki';  case 7    ede = 2; alfa = 4; model1 = 'PR-Yokozeki';endtemp = textscan(fid, '%s %d', 1,'HeaderLines', 8);poflag = cell2mat(temp(2));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 2);a0 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);a1 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);a2 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);a3 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);b0 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);b1 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);b2 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);b3 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);c0 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);c1 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);c2 =  cell2mat(temp(4));temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);
c3 =  cell2mat(temp(4));
temp = textscan(fid, '%s %s %s %f', 1,'HeaderLines', 0);
c4 =  cell2mat(temp(4));
temp = textscan(fid, '%s %d', 1,'HeaderLines', 7);moflag = cell2mat(temp(2));rflag = cell2mat(textscan(fid, '%d', 1,'HeaderLines', 2));gflag = cell2mat(textscan(fid, '%d', 1,'HeaderLines', 2));expdata = cell2mat(textscan(fid, '%f %f %f', 'HeaderLines', 2));fclose(fid);disp(' ')disp('Calculando ...') disp(' ')pexp = expdata(:, 1);pexp = 10 * pexp;texp = expdata(:, 2);etaexp = expdata(:, 3);np = size(expdata, 1);# calcular presiones de repulsion/atraccion (PR & PA) y densidad con la EdEC[pr, pa, rhocal] = prhoede(pexp, texp);  # armar vector con los estimados iniciales de los parametros a optimizarswitch poflag  case 1    x = [a1 a2 a3 b1 b2 b3 c2 c3]; model2 = 'TF8P-';  case 2    x = [a1 a2 a3 b1 b2 b3 c4]; model2 = 'TF7P-';  case 3    x = [a0 a1 a2 b0 b1 b2 c0 c1 c2]; model2 = 'TF9P-';  case 4    x = [a1 a2 b1 b2 c2 c3]; model2 = 'TF6P-';  case 5    x = [b1 c2 c3]; model2 = 'TF3P-';end# definir la funcion de minimos cuadrados (solo para el metodo de Levenberg-Marquardt)function y = flm(varind, x)   global pexp texp etaexp
  etacal = visco(texp, x);
  y = etacal ./ etaexp;endfunction# definir la funcion de minimos cuadradosfunction y = f(x)  global pexp texp etaexp  etacal = visco(texp, x);  y = sum((1 - etacal ./ etaexp).^2);endfunction# minimizacion de la funcion objetivo por varios metodosif rflag == 1  switch moflag    case 1                     # L-M      varobs = ones(size(etaexp));
      [fval, xfit, cvg, iter] = leasqr(texp, varobs, x, @flm);      xfit, cvg, iter    case 2                     # N-M       [xfit, fval, nev] = nelder_mead_min(@f, x)          case 3                     # POWELL       [xfit, fval, convergence, iters, nevs] = powell(@f, x)    case 4                     # D-E      switch poflag        case 1          ctl.XVmin = [-1 -.1 -.01 -1 -.1 -.01 -.1 -.01];          ctl.XVmax = [ 1  .1  .01  1  .1  .01  .1  .01];        case 2          ctl.XVmin = [-1 -.1 -.01 -1 -.1 -.01 -.01];
          ctl.XVmax = [ 1  .1  .01  1  .1  .01  .01];
        case 3          ctl.XVmin = [-1 -.1 -.01 -1 -.1 -.01 -.1 -.01 -.001];
          ctl.XVmax = [ 1  .1  .01  1  .1  .01  .1  .01  .001];
        case 4          ctl.XVmin = [-1 -.1 -.01 -1 -.1 -.01];          ctl.XVmax = [ 1  .1  .01  1  .1  .01];        case 5          ctl.XVmin = [-.01 -1e-5 -1e-8];          ctl.XVmax = [ .01  1e-5  1e-8];      end      [xfit, obj_value, nfeval, convergence] = de_min(@f, ctl)      end  disp(' ')else  xfit = x;endif# guardar resultados en el archivo de salida data_out.txtdisp('Guardando los resultados en data_out.txt ...')disp(' ')etacal = visco(texp, xfit);rdev_eta = 100 * (1 - etacal ./ etaexp);dap_eta = sum(abs(rdev_eta)) / np;bias_eta = sum(rdev_eta) / np;mda_eta = max(abs(rdev_eta));i = 1: np;resultados = [i; texp'; pexp'; rhocal'; etaexp'; etacal'; rdev_eta'];fid = fopen('data_out.txt', 'w');fprintf(fid, 'Optimizacion de los coeficientes de friccion en la teoria F en combinacion con una EdEC\n');fprintf(fid, 'para el modelado de viscosidades de LI & DES\n\n');fprintf(fid, 'Compuesto: %s\n\n', compuesto);fprintf(fid, '*** Ecuacion de estado cubica: ');switch ede  case 1    fprintf(fid, 'S-R-K\n');  case 2    fprintf(fid, 'P-R\n');endfprintf(fid, '*** Funcion Alfa: ');switch alfa  case 1    fprintf(fid, 'Soave\n');  case 2    fprintf(fid, 'Stryjek-Vera\n');  case 3    fprintf(fid, 'Melhem et al.\n');  case 4    fprintf(fid, 'Yokozeki\n');endif (rflag == 1)  fprintf(fid, '*** Metodo de optimizacion: ');  switch moflag    case 1      fprintf(fid, 'L-M\n');    case 2      fprintf(fid, 'N-M\n');    case 3      fprintf(fid, 'POWELL\n');    case 4      fprintf(fid, 'D-E\n');  endendiffprintf(fid, '\n');if (rflag == 1)  fprintf(fid, 'Parametro        Valor optimizado\n');  fprintf(fid, '---------------------------------\n');else  fprintf(fid, 'Parametro        Valor especificado\n');  fprintf(fid, '-----------------------------------\n');endifswitch poflag  case 1    fprintf(fid, 'A1 [cP/bar]       %14.7e\n', xfit(1));    fprintf(fid, 'A2 [cP/bar]       %14.7e\n', xfit(2));    fprintf(fid, 'A3 [cP/bar]       %14.7e\n', xfit(3));    fprintf(fid, 'B1 [cP/bar]       %14.7e\n', xfit(4));    fprintf(fid, 'B2 [cP/bar]       %14.7e\n', xfit(5));    fprintf(fid, 'B3 [cP/bar]       %14.7e\n', xfit(6));    fprintf(fid, 'C2 [cP/bar2]      %14.7e\n', xfit(7));    fprintf(fid, 'C3 [cP/bar2]      %14.7e\n\n\n', xfit(8));  case 2    fprintf(fid, 'A1 [cP/bar]       %14.7e\n', xfit(1));
    fprintf(fid, 'A2 [cP/bar]       %14.7e\n', xfit(2));
    fprintf(fid, 'A3 [cP/bar]       %14.7e\n', xfit(3));
    fprintf(fid, 'B1 [cP/bar]       %14.7e\n', xfit(4));
    fprintf(fid, 'B2 [cP/bar]       %14.7e\n', xfit(5));
    fprintf(fid, 'B3 [cP/bar]       %14.7e\n', xfit(6));
    fprintf(fid, 'C4 [cP/bar2]      %14.7e\n\n\n', xfit(7));
  case 3    fprintf(fid, 'A0 [cP/bar]       %14.7e\n', xfit(1));
    fprintf(fid, 'A1 [cP/bar]       %14.7e\n', xfit(2));
    fprintf(fid, 'A2 [cP/bar]       %14.7e\n', xfit(3));
    fprintf(fid, 'B0 [cP/bar]       %14.7e\n', xfit(4));
    fprintf(fid, 'B1 [cP/bar]       %14.7e\n', xfit(5));
    fprintf(fid, 'B2 [cP/bar]       %14.7e\n', xfit(6));
    fprintf(fid, 'C0 [cP/bar2]      %14.7e\n', xfit(7));
    fprintf(fid, 'C1 [cP/bar2]      %14.7e\n', xfit(8));
    fprintf(fid, 'C2 [cP/bar2]      %14.7e\n\n\n', xfit(9));
  case 4    fprintf(fid, 'A1 [cP/bar]       %14.7e\n', xfit(1));    fprintf(fid, 'A2 [cP/bar]       %14.7e\n', xfit(2));    fprintf(fid, 'B1 [cP/bar]       %14.7e\n', xfit(3));    fprintf(fid, 'B2 [cP/bar]       %14.7e\n', xfit(4));    fprintf(fid, 'C2 [cP/bar2]      %14.7e\n', xfit(5));    fprintf(fid, 'C3 [cP/bar2]      %14.7e\n\n\n', xfit(6));  case 5    fprintf(fid, 'B1 [cP/bar]       %14.7e\n', xfit(1));    fprintf(fid, 'C2 [cP/bar2]      %14.7e\n', xfit(2));    fprintf(fid, 'C3 [cP/bar2]      %14.7e\n\n\n', xfit(3));endfprintf(fid, '  N    Texp [K]   Pexp [bar]    RHOcal     ETAexp     ETAcal     %%ETAdev\n');fprintf(fid, '-------------------------------------------------------------------------\n');fprintf(fid, '%3d %10.2f %11.2f %11.2f %10.2f %10.2f %10.2f\n', resultados); fprintf(fid, '-------------------------------------------------------------------------\n');fprintf(fid, 'RHOcal en kg/m3\n');fprintf(fid, 'Todas las ETAs en mPa-s\n\n');fprintf(fid, '*** %%DAP-ETA  = %7.4f\n', dap_eta);fprintf(fid, '*** %%BIAS-ETA = %7.4f\n', bias_eta);fprintf(fid, '*** %%MDA-ETA  = %7.4f\n', mda_eta);fclose(fid);if (gflag == 1)# graficar los resultados obtenidos  disp('Imprimiendo graficas ...')  disp(' ')# grafica 3D de ETAexp & ETAcal vs T & P
  figure
  grid
  scatter3(pexp, texp, etaexp, 80, 'r', 'filled')
  xlabel('Pressure [bar]')
  ylabel('Temperature [K]')
  zlabel('Viscosity [mPa-s]')
  title(sprintf('Experimental & Calculated Viscosities vs T & P for %s', compuesto))
  hold on
  pp = linspace (min(pexp), max(pexp), 100)';  tt = linspace (min(texp), max(texp), 100)';  [xx, yy] = meshgrid(pp, tt);  [pr, pa, rhocal] = prhoede(xx, yy);  zz = visco(yy, xfit);
  surf(xx, yy, zz, 'FaceAlpha', .75, 'EdgeColor', 'none');
  set(gca,'ZScale','log')
  set(gcf, 'Position', [200, 200, 800, 800]);
  view(135, 25)
  set(gca, 'fontsize', 18);
  label1 = [model2, model1, ' Model'];
  legend('Experimental Data', label1)
  # grafica de RDEV_ETA vs P teniendo como parametro T
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
endifdisp('Ciao!')disp(' ')