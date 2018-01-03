% test_grafic.m
clc
clear all
 close all


%% Continuous case
x= linspace (-pi, pi, 1001);
f=exp(3*cos(x))/(2*pi*besseli(0,3));  
h=2*pi/1000;


     N=10;
    xp = linspace (-pi, pi, N + 1); % nodes
    fp=exp(3*cos(xp))/(2*pi*besseli(0,3)); % data 
    P=0; 
    

      % -----------------------------
       % case 1  % Lagrange
            for j=1:N+1;
                L=1;
                for i=1:N+1
                    if i~=j
                        L=conv(L,[1, -xp(i)]/(xp(j)-xp(i)));
                    end
                end
                P=P+fp(j)*L;
            end
           PN1=polyval(P,x);  
       % -----------------------------    
%         case 2 % Newton
            F=zeros(N+1,N+1);
            F(:,1)=fp;
            P=zeros(1,N+1);
            P(end)=F(1,1);
            
            for i=2:N+1
                P1=zeros(1,N-i+2); P1(end)=1;
                for j=2:i
                    F(i,j)=(F(i,j-1)-F(i-1,j-1))/(xp(i)-xp(i-j+1));
                    P1=conv(P1, [1, -xp(j-1)]);
                    
                end
                P=P+F(i,i)*P1;
            end
            PN2=polyval(P,x);
            
        %----------------------------------
%         case 3  %
            for i=1:N
                k=find((x>=xp(i))&(x<=xp(i+1)));
                PN3(k)=fp(i)+(fp(i+1)-fp(i))/(xp(i+1)-xp(i))*(x(k)-xp(i));
                
            end

    
%         case 4
            % --------------------------
            a=fp';
            h=(xp(2)-xp(1));
            
            A=sparse([1,N+1], [1,N+1],1,N+1,N+1)+sparse(2:N,1:N-1,h,N+1,N+1)...
                +sparse(2:N,3:N+1,h,N+1,N+1)+sparse(2:N, 2:N,4*h,N+1,N+1);
            B=zeros(N+1,1);
            i=2:N;  B(i)=3/h*(a(i+1)-2*a(i)+a(i-1));
            
            c=A\sparse(B);
            b=zeros(N+1,1); d=zeros(N+1,1);
            i=1:N;
            b(i)=(a(i+1)-a(i))/h-h*(c(i+1)+2*c(i))/3;
            d(i)=(c(i+1)-c(i))/(3*h);
            
            
            for j=1:N;
                k=find((x>=xp(j))&(x<=xp(j+1)));
                PN4(k)=a(j)+b(j)*(x(k)-xp(j))+c(j)*(x(k)-xp(j)).^2+d(j)*(x(k)-xp(j)).^3;
            end
            
            % ---------------------------------            
%         case 5

            a=fp';
            h=(xp(2)-xp(1));
            
            A=sparse([1,N+1], [1,N+1],2*h,N+1,N+1)+sparse(2:N+1,1:N,h,N+1,N+1)...
                +sparse(1:N,2:N+1,h,N+1,N+1)+sparse(2:N, 2:N,4*h,N+1,N+1);
            dfa=(fp(2)-fp(1))/h;dfb=(fp(N+1)-fp(N))/h;
            B=zeros(N+1,1);
            i=2:N;  B(i)=3/h*(a(i+1)-2*a(i)+a(i-1));
            B(1)=3/h*(a(2)-a(1))-3*dfa; B(N+1)=3*dfb-3/h*(a(N+1)-a(N));
            c=A\sparse(B);
            b=zeros(N+1,1); d=zeros(N+1,1);
            i=1:N;
            b(i)=(a(i+1)-a(i))/h-h*(c(i+1)+2*c(i))/3;
            d(i)=(c(i+1)-c(i))/(3*h);
            
            
            for j=1:N;
                k=find((x>=xp(j))&(x<=xp(j+1)));
                PN5(k)=a(j)+b(j)*(x(k)-xp(j))+c(j)*(x(k)-xp(j)).^2+d(j)*(x(k)-xp(j)).^3;
            end
  % ----------------------------------------------------           
%         case 6
            S=0;m=N/2;
            for k=1:m+1;
                a(k)=1/m*sum(fp.*cos(k*xp));
                b(k)=1/m*sum(fp.*sin(k*xp));
                
                if k==m+1;
                    S=S+a(k)*cos(k*x);
                else
                    S=S+a(k)*cos(k*x)+b(k)*sin(k*x);
                end
                    
            end
            a0=1/m*sum(fp);
            PN6=a0/2+S;
subplot(2,1,1);
plot(x,f,x,PN1,x,PN2,x,PN3,x,PN4,x,PN5,x,PN6)
title('Continous Case(N=10)')
xlabel('x'); ylabel('f(x)')
legend('f(x)','Lagrange','Newton', 'Linear', 'Natural','Cubic spline','Fourier')

%% Discreet Case

    Data=load('sunspot.dat');


    h=5;
    x =Data(:,1) ; % nodes
    f=Data(:,2); % sunspot
     h=h-1;

     I=1:h:300;
    xp=x(I); fp=f(I);
  
    N=length(xp)-1;
  
      % -----------------------------
%         case 1  % Lagrange
P=0;
            for j=1:N+1;
                L=1;
                for i=1:N+1
                    if i~=j
                        L=conv(L,[1, -xp(i)]/(xp(j)-xp(i)));
                    end
                end
                P=P+fp(j)*L;
            end
           PN1=polyval(P,x);  
           
       % -----------------------------    
%         case 2 % Newton
            F=zeros(N+1,N+1);
            F(:,1)=fp;
            P=zeros(1,N+1);
            P(end)=F(1,1);
            
            for i=2:N+1
                P1=zeros(1,N-i+2); P1(end)=1;
                for j=2:i
                    F(i,j)=(F(i,j-1)-F(i-1,j-1))/(xp(i)-xp(i-j+1));
                    P1=conv(P1, [1, -xp(j-1)]);
                    
                end
                P=P+F(i,i)*P1;
            end
            PN2=polyval(P,x);
            
        %----------------------------------
%         case 3  %
            for i=1:N
                k=find((x>=xp(i))&(x<=xp(i+1)));
                PN(k)=fp(i)+(fp(i+1)-fp(i))/(xp(i+1)-xp(i))*(x(k)-xp(i));
                
            end
            PN3=PN';
    
%         case 4
            % --------------------------
            a=fp;
            h=(xp(2)-xp(1));
            
            A=sparse([1,N+1], [1,N+1],1,N+1,N+1)+sparse(2:N,1:N-1,h,N+1,N+1)...
                +sparse(2:N,3:N+1,h,N+1,N+1)+sparse(2:N, 2:N,4*h,N+1,N+1);
            B=zeros(N+1,1);
            i=2:N;  B(i)=3/h*(a(i+1)-2*a(i)+a(i-1));
            
            c=A\sparse(B);
            b=zeros(N+1,1); d=zeros(N+1,1);
            i=1:N;
            b(i)=(a(i+1)-a(i))/h-h*(c(i+1)+2*c(i))/3;
            d(i)=(c(i+1)-c(i))/(3*h);
            
            
            for j=1:N;
                k=find((x>=xp(j))&(x<=xp(j+1)));
                PN(k)=a(j)+b(j)*(x(k)-xp(j))+c(j)*(x(k)-xp(j)).^2+d(j)*(x(k)-xp(j)).^3;
            end
            PN4=PN';
            % ---------------------------------            
%         case 5

            a=fp;
                       
            A=sparse([1,N+1], [1,N+1],2*h,N+1,N+1)+sparse(2:N+1,1:N,h,N+1,N+1)...
                +sparse(1:N,2:N+1,h,N+1,N+1)+sparse(2:N, 2:N,4*h,N+1,N+1);
            dfa=(fp(2)-fp(1))/h;dfb=(fp(N+1)-fp(N))/h;
            B=zeros(N+1,1);
            i=2:N;  B(i)=3/h*(a(i+1)-2*a(i)+a(i-1));
            B(1)=3/h*(a(2)-a(1))-3*dfa; B(N+1)=3*dfb-3/h*(a(N+1)-a(N));
            c=A\sparse(B);
            b=zeros(N+1,1); d=zeros(N+1,1);
            i=1:N;
            b(i)=(a(i+1)-a(i))/h-h*(c(i+1)+2*c(i))/3;
            d(i)=(c(i+1)-c(i))/(3*h);
            
            
            for j=1:N;
                k=find((x>=xp(j))&(x<=xp(j+1)));
                PN(k)=a(j)+b(j)*(x(k)-xp(j))+c(j)*(x(k)-xp(j)).^2+d(j)*(x(k)-xp(j)).^3;
            end
            PN5=PN';
  % ----------------------------------------------------           
%         case 6
            S=0;m=N/2;
            xp1=(xp-1850)/150*pi;
            x1=(x-1850)/150*pi;
            for k=1:m+1;
                a(k)=1/m*sum(fp.*cos(k*xp1));
                b(k)=1/m*sum(fp.*sin(k*xp1));
                
                if k==m+1;
                    S=S+a(k)*cos(k*x1);
                else
                    S=S+a(k)*cos(k*x1)+b(k)*sin(k*x1);
                end
                    
            end
            a0=1/m*sum(fp);
            PN6=a0/2+S;
x1=x(1:length(PN3));

subplot(2,1,2);
plot(x,f,x1,PN3,x1,PN4,x1,PN5,x,PN6)
title('Continous Case(N=60)')
xlabel('x'); ylabel('f(x)')
legend('f(x)', 'Linear', 'Natural','Cubic spline','Fourier')










