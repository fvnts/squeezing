%   ------------------------------------------------------------/
%   Ince–Strutt diagram of Mathieu equation
%
%   This is supplementary material for
%   https://arxiv.org/abs/1704.05840
%   coded by J Fuentes
%   jfuentes [at] fis.cinvestav.mx
%   ------------------------------------------------------------/
clc;
clear
close all;
%   ------------------------------------------------------------/
%   file
f = fopen('u12.dat','w');
g = fopen('c12.dat','w');
h = fopen('u21.dat','w');
k = fopen('c21.dat','w');
p = fopen('stab.dat','w');
q = fopen('nons.dat','w');
%   ------------------------------------------------------------/
%   interval
n = 256;
d = 2*pi/n;
%   ------------------------------------------------------------/
%   defs
J = [0,1;-1,0];
%   ------------------------------------------------------------/
%   iterative computation
for b0 = 0.50:1/n:3.5
    for b1 = -1.5:1/n:1.5
        
        U0 = eye(2);
        A  = zeros(2);
        
        for T = 0:d:2*pi
            b = b0 - 2*b1*sin(T); 
            A = [0,1;-b,0];
            U0= expm(d*A)*U0;
        end
            
        U = U0';
        
        %   trajectories
            
        if  (U(1,2) <= 9e-3 && U(1,2) >= 1e-21 )

            fprintf(f,'%f %f %f %f\n',U(1,:),U(2,:));
            fprintf(g,'%f %f\n',b0,b1);
            
        elseif (U(2,1) <= 9e-3 && U(2,1) >= 1e-21 )
            fprintf(h,'%f %f %f %f\n',U(1,:),U(2,:));
            fprintf(k,'%f %f\n',b0,b1);  
        end
        
        %   stability regions
        
        if      ( trace(U) <=2  && trace(U) >=-2 )
            fprintf(p,'%f %f\n',b0,b1);
        elseif  ( trace(U) > 2  || trace(U) < -2 )
            fprintf(q,'%f %f\n',b0,b1);
        end
        
    end
end

fclose(f); fclose(g); 
fclose(h); fclose(k);
fclose(p); fclose(q);

disp('done !')
%   ------------------------------------------------------------/
%   eof
%   ------------------------------------------------------------/
