%   ------------------------
%
%   squeezing probability
%
%   This is supplementary material for
%   https://arxiv.org/abs/1704.05840
%   coded by J Fuentes
%   jfuentes [at] fis.cinvestav.mx
%   ------------------------

clc; clear; close all;

%   ------------------------
%   open files

pr  = fopen('prob.dat','w');

dq  = importdata('shadow.dat');

q1  = importdata('trajectory1.dat'); 
q2  = importdata('trajectory2.dat');

%   :::::

dq1 = dq(1:200);
dq2 = flip( dq(201:end) ) * 0.875;

dq1 = dq1';
dq2 = dq2';

dq  = [dq1,dq2];  

%   :::::

q1  = q1(:,1);
q2  = q2(:,1);
q1  = q1';
q2  = q2';

q  = horzcat(q1,q2);

%   ------------------------
%   definitions

N =  256;
a = -2;
b =  4;
h = (b - a)/N;
x =  linspace(a,b,N+1);
r = zeros;
    
%   ------------------------
%   gaussian shadow

for k = 1:length(dq)

    m = 1;
    
    for j = 1:N + 1
    
        s =  sqrt(2) / ( sqrt(pi) * dq(k) ) * exp( - power( ( x(j) - q(k) ) ./ dq(k) , 2 )  );
    
        if s < 5e-1 && s > 1e-2
            
            r(m) = x(j);
            m    = m + 1;
            
        end
            
    end
    
    rmin = min(r)/2;
    rmax = max(r)/2;
    
    dr   = (rmax - rmin)/2;

%   CFL number >= 0.5 to be Liapunov–stable
    qmax = ( q(k) + dr*dq(k) ) + 0.45*dr;
    qmin = ( q(k) - dr*dq(k) ) - 0.45*dr; 
    
    if mod(k,7)
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
