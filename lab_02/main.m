% вариант 13

function main()
  close all;
  
  f = fopen("X.txt", "r");
  buf = textscan(f, "%f", 'Delimiter', ' ');
  X = buf{1,1};
  X = X';
  fclose(f);
  %X = sort(X);
  fprintf("%f\n", X);
  
  gamma = 0.9;
  alpha = (1 - gamma) / 2;
  
  % 1.1) и 2) вычисление точечных оценок MX и DX
  MU = MX(X);
  S2 = DX(X);
  fprintf("MU = %f, S2 = %f\n", MU, S2);
  
  % 1.2) границы доверительного интервала для MU
  MU_low = get_low_MU(X, alpha);
  MU_high = get_high_MU(X, alpha);
  fprintf("MU_low = %f, MU_high = %f\n", MU_low, MU_high);
  
  % 1.3) границы доверительного интервала для S2
  S2_low = get_low_S2(X, alpha);
  S2_high = get_high_S2(X, alpha);
  fprintf("S2_low = %f, S2_high = %f\n", S2_low, S2_high);
  
  % 3) пользовательский ввод и графики
  gamma = 0.9;
  alpha = (1 - gamma) / 2;
  N = length(X);
  start = 5;
  
  figure;
  hold on;
  plot([start, N], [MU, MU], 'r');
  MU_array = get_MU_array(X, N);
  plot((start : N), MU_array(start : N), 'g');
  MU_low_array = get_MU_low_array(X, N, alpha);
  plot((start : N), MU_low_array(start : N), 'b');
  MU_high_array = get_MU_high_array(X, N, alpha);
  plot((start : N), MU_high_array(start : N), 'm');
  l = legend('\mu\^(x_N)', '\mu\^(x_n)', '_{--}\mu\^(x_n)', '^{--}\mu\^(x_n)');
  set(l, 'fontsize', 14);
  grid on;
  
  figure;
  hold on;
  plot([start, N], [S2, S2], 'r');
  S2_array = get_S2_array(X, N);
  plot((start : N), S2_array(start : N), 'g');
  S2_low_array = get_S2_low_array(X, N, alpha);
  %fprintf("%f\n", S2_low_array);
  plot((start : N), S2_low_array(start : N), 'b');
  S2_high_array = get_S2_high_array(X, N, alpha);
  %fprintf("%f\n", S2_high_array);
  plot((start : N), S2_high_array(start : N), 'm');
  l = legend('S^2(x_N)', 'S^2(x_n)', '_{--}\sigma^2(x_n)', '^{--}\sigma^2(x_n)');
  set(l, 'fontsize', 14);
  grid on;
 end
 
 function mx = MX(X)
   n = length(X);
   sum = sum(X);
   mx = sum / n;
 endfunction
 
 function dx = DX(X)
   n = length(X);
   MX = MX(X);
   if (n == 1)
     dx = 0;
   else
     dx = sum((X - MX).^2) / (n - 1);
   endif
 endfunction
 
 function mu_low = get_low_MU(X, alpha)
   n = length(X);
   MU = MX(X);
   S = sqrt(DX(X));
   mu_low = MU + S / sqrt(n) * tinv(alpha, n - 1);
 endfunction
 
  function mu_high = get_high_MU(X, alpha)
   n = length(X);
   MU = MX(X);
   S = sqrt(DX(X));
   mu_high = MU + S / sqrt(n) * tinv(1 - alpha, n - 1);
 endfunction
 
 function s2_low = get_low_S2(X, alpha)
   n = length(X);
   S2 = DX(X);
   s2_low = (n - 1) * S2 / chi2inv(1 - alpha, n -1);
 endfunction
 
 function s2_high = get_high_S2(X, alpha)
   n = length(X);
   S2 = DX(X);
   s2_high = (n - 1) * S2 / chi2inv(alpha, n -1);
 endfunction
 
 function MU_array = get_MU_array(X, N)
   MU_array = zeros(1, N);
   for i = 1 : N
     MU_array(i) = MX(X(1 : i));
   endfor
 endfunction
 
 function MU_low_array = get_MU_low_array(X, N, alpha)
   MU_low_array = zeros(1, N);
   for i = 1 : N
     cur_X = X(1 : i);
     MU_low_array(i) = get_low_MU(cur_X, alpha); 
   endfor
 endfunction
 
 function MU_high_array = get_MU_high_array(X, N, alpha)
   MU_high_array = zeros(1, N);
   for i = 1 : N
     cur_X = X(1 : i);
     MU_high_array(i) = get_high_MU(cur_X, alpha); 
   endfor
 endfunction
 
 function S2_array = get_S2_array(X, N)
   S2_array = zeros(1, N);
   for i = 1 : N
     S2_array(i) = DX(X(1 : i));
   endfor
 endfunction
 
 function S2_low_array = get_S2_low_array(X, N, alpha)
   S2_low_array = zeros(1, N);
   for i = 1 : N
     cur_X = X(1 : i);
     S2_low_array(i) = get_low_S2(cur_X, alpha);
   endfor
 endfunction
 
 function S2_high_array = get_S2_high_array(X, N, alpha)
   S2_high_array = zeros(1, N);
   for i = 1 : N
     cur_X = X(1 : i);
     S2_high_array(i) = get_high_S2(cur_X, alpha);
   endfor
 endfunction
