%   ------------------------
%
%   squeezing probability
%
%   This is supplementary material for
%   https://arxiv.org/abs/1704.05840
%   coded by J Fuentes
%   jfuentes [at] fis.cinvestav.mx
%
%   ------------------------

clc; clear; close all;

%   ------------------------
%   open files

pr  = fopen('prob.dat','w');
dq  = importdata('shadow.dat');
q1  = importdata('trajectory1.dat'); 

%   :::::

q1  = q1(:,1);
q1  = q1';

q  = q1;

%   ------------------------
%   definitions

N = 800;
a = -2;
b = 4;
h = (b - a)/N;
x = linspace(a,b,N+1);
r = zeros;
    
%   ------------------------
%   gaussian shadow

for k = 1:length(dq)

    m = 1;
    
    for j = 1:N + 1
    
        s =  sqrt(2) / ( sqrt(pi) * dq(k) ) * exp( - power( ( x(j) - q(k) ) ./ dq(k) , 2 )  );
    
        if s < 6e-1 && s > 1e-5
            
            r(m) = x(j);
            m    = m + 1;
            
        end
            
    end
    
    rmin = min(r)/2;
    rmax = max(r)/2;
    
    dr   = (rmax - rmin)/2;

    %   numerical tolerance is set >= 0.15 to be Liapunov–stable
    qmax = ( q(k) + dr*dq(k) ) + 0.15*dr;
    qmin = ( q(k) - dr*dq(k) ) - 0.15*dr;
    
    if mod(k,5)
       fprintf(pr,'%d %f %f \n', k, NaN, NaN);  
    else
       fprintf(pr,'%d %f %f \n', k, qmin, qmax); 
    end
end

%   ------------------------
%   test-plots

%plot(q)

%   ------------------------

fclose(pr);

%   eof
