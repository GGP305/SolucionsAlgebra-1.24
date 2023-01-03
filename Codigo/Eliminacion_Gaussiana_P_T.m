function [ret] = Eliminacion_Gaussiana_P_T(A,b)

							%Muestra una  respuesta con mayor ecxatitud y muestra la matriz aumentada
[n,m]=size(A);							% TamaÃ±o de la matriz A
Ab=[A b];								% Matriz aumentada [Ab].

if n==m 								% si la matriz es cuadrada
    for i=1:n
        marca(i)=i; 							%Almacena posiciones de las variables x en cada cambio
            
    end
    for k=1:(n-1)
        mayor=0; 							%Hace el mayor de la fila igual a cero
        filam=k; 							%La fila k contiene el numero mayor
        columnam=k; 							%Hace la columna k como la que tiene el numero mayor de toda la matriz
        for p=k:n
            for r=k:n
                if mayor<abs(Ab(p,r))				%Se busca el numero mayor por toda la matriz
                mayor=abs(Ab(p,r)); 					%Hace como el mayor el numero encontrado
                filam=p;							%cambio de fila del numero mayor
                columnam=r;						%cambio de columna del numero mayor
                end
            end
        end
         if mayor ==0
            fprintf(' El sistema tiene infinitas soluciones ')
            break							%se interrumpe el programa con la instruccion break, ya que 
                 								%si mayor=0, mas adelante se obtiene una division por cero.       								 
         else
           if filam ~= k                     %si el mayor no esta en la fila k
            for j=1:(n+1)
                aux=Ab(k,j);						%Se utiliza una variable auxiliar para intercambio de fila
                           							 
                Ab(k,j)=Ab(filam,j);
                Ab(filam,j)=aux;
            end
           end
            if columnam ~= k                  %si el mayor no esta en la columna k
            for i=1:n
                aux=Ab(i,k);						%Se utiliza una variable auxiliar para intercambio de columna				
                Ab(i,k)=Ab(i,columnam);
                Ab(i,columnam)=aux;
            end
            aux = marca(k);						%Se utiliza una variable auxiliar para intercambio de variables
            marca(k)= marca(columnam);
            marca(columnam)=aux;
            end
         end
 
         
         for i=(k+1):n
            mult(i,k)=Ab(i,k)/Ab(k,k); 					%Multiplicadores de cada etapa

            for j=k:(n+1)
                Ab(i,j)= Ab(i,j) - mult(i,k)*Ab(k,j);			 %Formula para convertir fila
            end
         end
    end
            marca(columnam)=aux;
     %Sustitucion regresiva
	 for i=n:-1:1
            suma=0;
               for p=(i+1):n
                suma = suma + Ab(i,p)*X(p);
               end
            X(i)=(Ab(i,n+1)-suma)/Ab(i,i);				
             
     end
     %la siguiente parte del programa ordena las varibles, tomando en
     %cuenta la marca final y los retoma con su coeficiente a la marca
     %inicial
         for i=1:n
             for j=1:n
                 if marca(j)==i
                     k=j;
                 end
             end
             aux=X(k);							%para poder intercambiar las variables, se utiliza una variable auxiliar.                   
             X(k)=X(i);
             X(i)=aux;
             aux=marca(k);						%para poder intercambiar las variables, se utiliza una variable auxiliar.
             marca(k)=marca(i);
             marca(i)=aux;
         end
else									 %En caso de que la matriz no sea cuadrada
     fprintf('No se puede realizar debido a que la matriz no es cuadrada ');
end

%Muestra los resultados
ret = transpose(X);