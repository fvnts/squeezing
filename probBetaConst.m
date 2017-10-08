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
q2  = importdata('trajectory2.dat');

%   :::::

dq1 = dq(1:200);
dq2 = flip( dq(201:400) );
dq3 = flip( dq(401:end) )*1.5;

dq1 = dq1';
dq2 = dq2';
dq3 = dq3';

dq  = [dq1,dq2,dq3];  

%   :::::

q1  = q1(:,1);
q2  = q2(:,1);
q1  = q1';
q2  = q2';

q  = horzcat(q1,q2);

%   ------------------------
%   definitions

b0 = 1/10;

N =  256;
a = -2;
b =  4;
h = (b - a)/N;
x =  linspace(a,b,N+1);
t =  linspace(-pi/2,pi/(2*sqrt(b0)),length(q));
r = zeros;
    
%   ------------------------
%   gaussian shadow

for k = 1:length(dq)

    m = 1;
    
    for j = 1:N + 1
    
        s =  sqrt(2) / ( sqrt(pi) * dq(k) ) * exp( - power( ( x(j) - q(k) ) ./ dq(k) , 2 )  );
    
        if s < 5e-1 && s > 5e-3
            
            r(m) = x(j);
            m    = m + 1;
            
        end
            
    end
    
    rmin = min(r);
    rmax = max(r);
    
    dr   = (rmax - rmin)/2.5;

%   code is stable in the whole space
    qmax = q(k) + dr;
    qmin = q(k) - dr; 
    
    if mod(k,7)
       fprintf(pr,'%d %f %f \n', k, NaN, NaN);  
    else
       fprintf(pr,'%d %f %f \n', k, qmin, qmax); 
    end
end

%   ------------------------
%   test-plots

%plot(t,q)

%   ------------------------

fclose(pr);

%   eof
