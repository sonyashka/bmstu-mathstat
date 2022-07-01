% вариант 13

function main()
  f = fopen("X.txt", "r");
  buf = textscan(f, "%f", 'Delimiter', ' ');
  X = buf{1,1};
  X = X';
  fclose(f);
  %X = [-10.82,-9.27,-9.65,-9.36,-9.27,-11.25,-9.89,-9.26,-11.15,-8.90,-11.02,-8.28,-9.18,-10.16,-10.59,-10.82,-9.05,-9.47,-10.98,-11.50,-8.64,-10.86,-10.76,-11.49,-11.09,-9.33,-9.32,-9.66,-8.79,-10.54,-9.12,-10.40,-8.59,-10.22,-9.06,-10.59,-10.60,-10.25,-9.35,-11.44,-9.81,-9.32,-9.95,-9.33,-10.64,-9.45,-10.99,-10.15,-10.39,-10.36,-10.49,-11.67,-10.00,-10.87,-11.11,-9.68,-10.77,-9.13,-8.62,-10.33,-11.36,-10.24,-9.41,-11.05,-10.15,-9.35,-11.45,-9.87,-10.41,-10.11,-10.84,-11.48,-7.77,-10.79,-9.88,-10.70,-9.07,-9.47,-10.15,-9.93,-11.52,-9.04,-10.93,-10.13,-9.56,-11.39,-9.79,-9.19,-11.09,-9.86,-10.67,-10.26,-9.07,-10.53,-11.24,-10.16,-11.33,-8.76,-8.88,-10.53,-10.12,-8.98,-9.84,-9.90,-10.13,-9.32,-9.31,-9.99,-8.55,-11.64,-11.32,-10.51,-11.71,-10.50,-10.50,-12.20,-11.68,-10.45,-7.88,-10.84];
  X = sort(X);
  fprintf("%f\n", X);

  % а) M_max и M_min
  M_min = X(1);
  M_max = X(end);
  fprintf("M_min = %d, M_max = %d\n", M_min, M_max);

  % б) R
  R = M_max - M_min;
  fprintf("R = %d\n", R);

  % в) исправленная дисперсия
  MU = MX(X);
  S2 = DX(X);
  fprintf('MU = %d, S^2 = %d\n', MU, S2);

  % г) группирорвка значений выборки в m интервала
  m = floor(log2(length(X)) + 2);
  groupOnIntervals(X, m);
  hold on;

  sigma = sqrt(S2);
  % д) диаграмма и функция плотности распределения
  x_values = (M_min - 1.5):(R / 100):(M_max + 1.5);
  Y = zeros(1, length(x_values));
  for i = 1:length(x_values)
    Y(i) = (1 / (sqrt(2 * pi) * sigma)) * exp(-(x_values(i) - MU) * (x_values(i) - MU) / (2 * sigma * sigma));
  endfor
  plot(x_values, Y), grid;
  
  % е) эмпирическая функция распределения и функция распределения (норм.)
  figure;
  
  Y = zeros(1, length(X)+10);
  Xd = zeros(1, length(X)+10);
  for i = 1 : 5
    Xd(i) = X(1) - (6-i)*0.1;
    Xd(i+length(X)+5) = X(end) + 0.1 * i;
  endfor
  for i = 1 : length(X)
    Xd(i+5) = X(i);
  endfor
  
  for i = 1 : length(Xd)
    if Xd(i) < X(1)
      Y(i) = 0;
    elseif Xd(i) > X(length(X))
      Y(i) = 1;
    else
      Y(i) = (i - 5) / length(X);
    endif
  endfor
  stairs(Xd, Y), grid;
  
  hold on;
  
  Y = zeros(1, length(x_values));
  func = @(x) exp(-(x - MU).* (x - MU) / (2 * S2));
  for i = 1:length(x_values)
    Y(i) = 1 / (sigma * sqrt(2 * pi)) * integral(func, -Inf, x_values(i));
  endfor
  plot(x_values, Y), grid;
 end
 
 function mx = MX(X)
   n = length(X);
   sum = sum(X);
   mx = sum / n;
 endfunction
 
 function dx = DX(X)
   n = length(X);
   MX = MX(X);
   dx = sum((X - MX).^2) / (n - 1);
 endfunction
 
 function groupOnIntervals(X, m)
   n = length(X);
   delta = (X(end) - X(1)) / m;
   intervals = zeros(2, m + 1);
   intervals(1, 1) = X(1);
   for i = 1:m
     intervals(1, i + 1) = intervals(1, i) + delta; 
   endfor
   
   interval_ind = 1;
   for i = 1:n
     while (X(i) >= intervals(1, interval_ind + 1))
       interval_ind += 1;
     endwhile
     intervals(2, interval_ind) += 1;
   endfor
   
   fprintf("Кол-во интервалов разбиения: %d\n", m);
   for i = 1:m-1
     fprintf("[%6.2f; %6.2f) - %d значений\n", intervals(1,i), intervals(1, i + 1), intervals(2, i));
   endfor
   fprintf("[%6.2f; %6.2f] - %d значений\n", intervals(1, m), intervals(1, m + 1), intervals(2, m));
   
   for i = 1:m
     intervals(2, i) /= (n * delta);
   endfor
   
   stairs(intervals(1, :), intervals(2, :)), grid;
 endfunction
