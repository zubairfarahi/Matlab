function N=eval_interpolator_c(tip, eps);
% Eval_interpolator - determine how quickly a polynomial converges Interpolation  
% Tip -   1 lagrange  
%            2 newton 
%            3 linear spline  
%            4 natural  
%            5 cubic spline 
%            6 fourrier  
% Eps - accepted tolerance for convergence


x= linspace (-pi, pi, 1001);
f=exp(3*cos(x))/(2*pi*besseli(0,3));  
h=2*pi/1000;
E=1; N=10; dE=1;


while (E>eps)
     N=N+2;
    xp = linspace (-pi, pi, N + 1); % nodes
    fp=exp(3*cos(xp))/(2*pi*besseli(0,3)); % data 
    P=0; 
    
    switch tip
      % -----------------------------
        case 1  % Lagrange
            for j=1:N+1;
                L=1;
                for i=1:N+1
                    if i~=j
                        L=conv(L,[1, -xp(i)]/(xp(j)-xp(i)));
                    end
                end
                P=P+fp(j)*L;
            end
           PN=polyval(P,x);  
       % -----------------------------    
        case 2 % Newton
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
            PN=polyval(P,x);
            
        %----------------------------------
        case 3  %
            for i=1:N
                k=find((x>=xp(i))&(x<=xp(i+1)));
                PN(k)=fp(i)+(fp(i+1)-fp(i))/(xp(i+1)-xp(i))*(x(k)-xp(i));
                
            end

    
        case 4
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
                PN(k)=a(j)+b(j)*(x(k)-xp(j))+c(j)*(x(k)-xp(j)).^2+d(j)*(x(k)-xp(j)).^3;
            end
            
            % ---------------------------------            
        case 5

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
                PN(k)=a(j)+b(j)*(x(k)-xp(j))+c(j)*(x(k)-xp(j)).^2+d(j)*(x(k)-xp(j)).^3;
            end
  % ----------------------------------------------------           
        case 6
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
            PN=a0/2+S;
    end

E=[E, (h*sum(abs(f-PN).^2))^(1/2)];

if length(E)>10
    if E(end)-E(end-1)>0
        disp('The interrogator does not converge.')
        N=inf;
        return;
    end
    dE=abs(E(end)-E(end-1));
end

end

