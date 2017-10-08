%   ------------------------------------------------------------/
%   uncertainty shadows
%
%   This is supplementary material for
%   https://arxiv.org/abs/1704.05840
%   coded by J Fuentes
%   jfuentes [at] fis.cinvestav.mx
%   ------------------------------------------------------------/

clc; clear; close;

%   --------------

beta0 = 1/10;

%   files

pulse1      = fopen('beta1.dat','w');
pulse2      = fopen('beta2.dat','w');

trajectory1 = fopen('trajectory1.dat','w');
trajectory2 = fopen('trajectory2.dat','w');

sha = fopen('shadow.dat','w');
ash = fopen('antiShadow.dat','w');

heisenberg  = fopen('heisenberg.dat','w');

%   interval

s   = 200;
d   = pi/s;
m   = 1;
k   = 1;

%   ini

q = 1;
p = 1;

x0 = [q;p]; 
x1 = [0;0];
x2 = [0,0];

l  = length(-pi/2:d:pi/2-d);

shadow      = zeros(1,2*l);
antishadow  = zeros(1,2*l);

%   --------------------------------
%   beta pulse 1
%   --------------------------------

U0 = eye(2);
A  = zeros(2);
n  = 0;

%   --------------------------------

for t = -pi/2:d:pi/2-d

    n  = n + 1;
    
    U  = U0;
    x1 = U*x0;
    
    shadow(n)      = sqrt(1/(2*k)*U(1,1)^2 + (k/2)*U(1,2)^2);
    antishadow(n)  = sqrt((k/2)*U(2,2)^2 + 1/(2*k)*U(2,1)^2);

    fprintf(sha,'%f \n',shadow(n));
    fprintf(ash,'%f \n',antishadow(n));
    
    fprintf(pulse1,'%d %f %f %f %f\n',n,U(1,:),U(2,:));
    fprintf(trajectory1,'%f %f\n',x1(1),x1(2));
    
    n1 = 6654979829 + 3254716032*cos(2*t) - 1940023077*cos(4*t) ...
                 - 952388778*cos(6*t) - 11277513*cos(8*t) + 2849706*cos(10*t) ...
                 + 423801*cos(12*t);
    d1 = 8*(57349 + 2313*cos(2*t) - 369*cos(4*t) - 93*cos(6*t))^2;
          
    beta1 = n1./d1;
    
    A  = [0,1;-beta1,0];
    U0 = expm(d*A)*U;

end

%   --------------------------------
%   beta pulse 2
%   --------------------------------

U0  = eye(2);
A   = zeros(2);

n   = 0;
k   = sqrt(beta0);
        
%   --------------------------------

for t = 0:d:pi/(2*k)-d
            
    n  = n + 1;
    
    U  = U0;
  	x2 = U*x1;
    
    shadow(n+l)      = sqrt(1/(2*k)*U(1,1)^2 + (k/2)*U(1,2)^2);
    antishadow(n+l)  = sqrt((k/2)*U(2,2)^2 + 1/(2*k)*U(2,1)^2);
    
    fprintf(sha,'%f \n',shadow(n+l));
    fprintf(ash,'%f \n',antishadow(n+l));

    fprintf(pulse2,'%d %f %f %f %f\n',n,U(1,:),U(2,:));
    fprintf(trajectory2,'%f %f\n',x2(1),x2(2));
    
    beta2 = beta0;
    
    A   =  [0,1;-beta2,0];
    U0  =  expm(d*A)*U;
        
end

%   --------------------------------
         
for n = 1:2*l            
    uncert  = shadow(n)*antishadow(n);
    fprintf(heisenberg,'%f \n',uncert);
end

%   --------------------------------
        
fclose(pulse1);
fclose(pulse2);

fclose(trajectory1);
fclose(trajectory2);

fclose(sha);
fclose(ash);

fclose(heisenberg);

disp('done !')

%   eof
