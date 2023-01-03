function [x_grid, t_grid, U] = solve_crank_nicholson_heat_equation(L, T, I, J)
  %% Resuelve la ecuacion de calor en derivadas parciales usando la discretizacion 
  %% del dominio y las diferencias finitas. Para cada punto (x_i, t_i) del
  %% dominio, se genera el elemento de la matriz U en la posicion i y j , que 
  %% contiene el valor de U(x_i, t_j).
  %% 
  %% Input:
  %%    - L: Longitud del intervalo de definicion de la variable espacial.
  %%    - T: Longitud del intervalo de definición de la variable temporal.
  %%    - I: Numero de particiones que se realizan sobre el intervalo [0, L].
  %%    - J: Numero de particiones que se realizan sobre el intervalo [0, T].
  %% 
  %% Output:
  %%    - x_grid: Un vector, representando la particion efectuada en el eje X.
  %%    - t_grid: Un vector, representando la particion efectuada en el eje Y.
  %%    - U: Una matriz, representando los valores de la funcion discretizada.
  
  %%%%%%%%%%%%%%%%%% Parametros iniciales %%%%%%%%%%%%%%%%%%%%
  x_grid = linspace(0, L, I+1);
  t_grid = linspace(0, T, J+1);
  U = zeros(I+1, J+1);
  
  h = L/I
  k = T/J
  alpha = k/h**2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%%%%%%%%%% Inicializacion de la matriz del sistema %%%%%%%%
  B = zeros(I-1, I-1);
  
  for k=1:I-1
    % Diagonal principal
    B(k,k) = 1 + alpha;
  endfor
  B
  for k=2:I-1
    % Subdiagonal superior
    B(k-1,k) = -alpha/2;
  endfor
  B
  for k=1:I-2
    % Subdiagonal inferior
    B(k+1,k) = -alpha/2;
  endfor
  B
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%% Condiciones de frontera %%%%%%%%%%%%%%%%%%
  for i = 1:I+1
    U(i, 1) = x_grid(i)*(x_grid(i)-1);
  endfor
  U
  for j = 1:J+1
    U(1, j) = t_grid(j);
    U(I+1, j) = -t_grid(j);
  endfor
  U
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%% Resolucion iterativa del sistema de ecuaciones %%%%%%
  for k = 2:J+1
    independent_term = 2*U(2:I, k-1);
    independent_term(1) = 2*U(1, k-1) + (alpha*(U(1,k)+U(1,k-1)))/2;
    independent_term(I-1) = 2*U(I, k-1) + (alpha*(U(I+1,k)+U(I+1,k-1)))/2;
    
    [Ws] = eliminacion_gaussiana_p_t(B,independent_term);
    Ws
    U(2:I, k) = Ws-U(2:I,k-1);
  endfor
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  U
endfunction
