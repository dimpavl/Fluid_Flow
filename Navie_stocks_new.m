clc;
clear;
%����� ����������
img = imread('Card.bmp');
 
%����������� ���������
V = 100; % ��������
p = 10^3; % ���������
M0 = 10^(-3); %��������
R = 100*10^(-9);%������� ����� 
D = 10^(-6);%����������� ��������
nu = M0/p; %�������������� ��������
Re = R*V/nu;%����� ����������
Prd = R*V/D;%�������� ������������
P = 10^6; % ����������� �������� 1 ���������
T = R/V; %������� �� �������
%��������� ��� ��������� ������������ ��������
alpha = 1;

Ce = 0.1;
C0 = 1;



m0 = 1;
u = 1;% ����������� ������� ������
P0 = 0.1;%�������� �� ����� �������������
gamma = 10^(floor(log10(P0)));

%��������� �����
dx =10^(-2); 
H = 1; L = 0.4;
A = 10^(3);

% ��������� �������
timeend = 7.1;


N = floor(H/dx);%���������� ������
M = floor(L/dx);%���������� ������
Vxold(M, N) = 0;
Vyold(M, N) = 0;
Pold(M, N) = 0;
Vxnew(M, N) = 0;
Vynew(M, N) = 0;
Pnew(M, N) = 0; 
Cnew(M, N) = 0;
Cold(M, N) = 0;
m(M, N) = 0;
DeltaP(M, N) = 0;

%����� ��������� ��� N
Nold(M, N) = 0;
Nnew(M, N) = 0;
K = 1/Prd;
Nx = 1/(dx^2);
Hr = 1;
%��� ���������� ������� �� ��������
l = 0

[X,Y] = MeshGrid(1:1:N,1:1:M);

fprintf(1,'��������� ��������\n');
fprintf(1,'-----------------------------------------\n');
fprintf(1,'��� ������� �������� �� ������������\n');
fprintf(1,'-----------------------------------------\n');
fprintf(1,'alpha\n');
alpha
fprintf(1,'Prd �������� ������������\n');
Prd
fprintf(1,'Ce\n');
Ce
fprintf(1,'��������� ������������ C0\n');
C0
fprintf(1,'-----------------------------------------\n');
fprintf(1,'��� ������� ��������� �����-������\n');
fprintf(1,'-----------------------------------------\n');
fprintf(1,'����������� �������������� ��������(���������)\n');
m0

fprintf(1,'��� �����\n');
dx
fprintf(1,'�������� �� ����� ������\n');
P0
fprintf(1,'���������� ������\n');
H
fprintf(1,'���������� ������\n');
L
fprintf(1,'gamma \n');
gamma
fprintf(1,'-----------------------------------------\n');
fprintf(1,'��������� ����������\n');
[0 timeend]
fprintf(1,'-----------------------------------------\n');
if (u == 0)
    fprintf(1,'���������� ��������\n');
else
    fprintf(1,'���������� ��������\n');
end;


%��������� �������
   for i = 1:1:M     
      if (granitsy(i,1,dx,L,H,img) == 1) 
        Vynew(i,1) = A*((i-1)*dx)^(2)*(L-(i)*dx)^(2);
        Vxnew(i,1)=0;
        Pnew(i,1) = P0;
        Cnew(i,1) = C0;        
      end;
      if (granitsy(i,1,dx,L,H,img) == 1)
        Pnew(i,N) = 0;       
      end;      
   end;
    
  for i = 1:1:M 
        for j = 1:1:N
            m(i,j) = m0;            
        end;
  end;
  
  
time = 0;
timek = 0;
k = 0;
flagexit = 1;
while (flagexit == 1)   
    if (time > timeend) flagexit = 0; end;
    if (flagexit == 0) continue; end;
    
    maxM = -Inf; 
    maxV = -Inf; 
    maxC = -Inf; 
    maxP = -Inf;
    maxN = -Inf;
    for i = 1:1:M 
        for j = 1:1:N 
           if(granitsy(i,j,dx,L,H,img) == 1)
               if (time <= timeend)
                    m(i,j) = m0*(1+u*Cold(i,j));                
                    t = sqrt(Vxold(i,j)^2+Vyold(i,j)^2);                
                    if (t > maxV) maxV = t; end;
                    if (m(i,j)>maxM) maxM = m(i,j); end;
                    if (Cold(i,j) > maxC) maxC = Cold(i,j); end;
                    if (Pold(i,j) > maxP) maxP = Pold(i,j); end;
                    if (Nold(i,j)> maxN) maxN = Nold(i,j); end;
               end;               
           end;
        end;
    end;
    
    if (time <= timeend)
        dt=0.1*min([dx/maxV, dx*dx/(nu)*0.1, dx^2/(1/Prd*maxC), dx^2/(K*maxN)]);
    end;    
    
    dt
    
    for i = 2:1:M-1 
        for j = 2:1:N-1  
            if (granitsy(i,j,dx,L,H,img) == 1)
                if (time <= timeend)                    
                    Vxnew(i,j) =  dt*(-0.5*(Vxold(i,j)+abs(Vxold(i,j))*(Vxold(i,j)-Vxold(i-1,j))/dx...
                                -0.5*(Vxold(i,j)-abs(Vxold(i,j)))*(Vxold(i+1,j)-Vxold(i,j))/dx...
                                -0.5*(Vyold(i,j)+abs(Vyold(i,j)))*(Vxold(i,j)-Vxold(i,j-1))/dx...
                                -0.5*(Vyold(i,j)-abs(Vyold(i,j)))*(Vxold(i,j+1)-Vxold(i,j))/dx))...             
                                +1/Re*m(i,j)*dt*(Vxold(i+1,j)+Vxold(i-1,j)+Vxold(i,j+1)+Vxold(i,j-1)-4*Vxold(i,j))/(dx)^2....
                                -1/Re*dt*P/(M0*V/R)*(Pold(i+1,j)-Pold(i-1,j))/(2*dx)...
                                +1/Re*2*dt*(m(i+1,j)-m(i-1,j))/(2*dx)*(Vxold(i+1,j)-Vxold(i-1,j))/(2*dx)...
                                +1/Re*dt*(m(i,j+1)-m(i,j-1))/(2*dx)*(Vyold(i+1,j)-Vyold(i-1,j)+Vxold(i,j+1)-Vxold(i,j-1))/(2*dx)...    
                                +Vxold(i,j);
  
                     Vynew(i,j) = dt*(-0.5*(Vxold(i,j)+abs(Vxold(i,j)))*(Vyold(i,j)-Vyold(i-1,j))/dx...
                                -0.5*(Vxold(i,j)-abs(Vxold(i,j)))*(Vyold(i+1,j)-Vyold(i,j))/dx...
                                -0.5*(Vyold(i,j)+abs(Vyold(i,j)))*(Vyold(i,j)-Vyold(i,j-1))/dx...
                                -0.5*(Vyold(i,j)-abs(Vyold(i,j)))*(Vyold(i,j+1)-Vyold(i,j))/dx)...
                                +1/Re*m(i,j)/p*dt/(dx)^2*(Vyold(i+1,j)+Vyold(i-1,j)+Vyold(i,j+1)+Vyold(i,j-1)-4*Vyold(i,j))...
                                -1/Re*dt*P/(M0*V/R)*(Pold(i,j+1)-Pold(i,j-1))/(2*dx)...
                                +1/Re*dt*(m(i+1,j)-m(i-1,j))/(2*dx)*(Vyold(i+1,j)-Vyold(i-1,j)+Vxold(i,j+1)-Vxold(i,j-1))/(2*dx)...
                                +1/Re*2*dt*(m(i,j+1)-m(i,j-1))/(2*dx)*(Vyold(i,j+1)-Vyold(i,j-1))/(2*dx)...
                                +Vyold(i,j);
                
                    Pnew(i,j) = -gamma*dt/(2*dx)*(Vxold(i+1,j)-Vxold(i-1,j)+Vyold(i,j+1)-Vyold(i,j-1))+Pold(i,j);
                
                    Cnew(i,j) = dt*1/Prd*(Cold(i+1,j)+Cold(i-1,j)+Cold(i,j+1)+Cold(i,j-1)-4*Cold(i,j))/(dx^2)...
                                +alpha*dt*(Ce-Cold(i,j))+Cold(i,j)+dt*...
                                (-0.5*(Vxold(i,j)+abs(Vxold(i,j))*(Cold(i,j)-Cold(i-1,j))/dx...
                                -0.5*(Vxold(i,j)-abs(Vxold(i,j)))*(Cold(i+1,j)-Cold(i,j))/dx...
                                -0.5*(Vyold(i,j)+abs(Vyold(i,j)))*(Cold(i,j)-Cold(i,j-1))/dx...
                                -0.5*(Vyold(i,j)-abs(Vyold(i,j)))*(Cold(i,j+1)-Cold(i,j))/dx));
                    
                    Nnew(i,j) = K*dt*(1-Nold(i,j)/Nx)*Cold(i,j)-dt*Hr*Nold(i,j)+Nold(i,j)+dt*...
                                (-0.5*(Vxold(i,j)+abs(Vxold(i,j))*(Nold(i,j)-Nold(i-1,j))/dx...
                                -0.5*(Vxold(i,j)-abs(Vxold(i,j)))*(Nold(i+1,j)-Nold(i,j))/dx...
                                -0.5*(Vyold(i,j)+abs(Vyold(i,j)))*(Nold(i,j)-Nold(i,j-1))/dx...
                                -0.5*(Vyold(i,j)-abs(Vyold(i,j)))*(Nold(i,j+1)-Nold(i,j))/dx));        
                end;                
            end;
        end;
    end;
    
%��������� �������
    for c = 1:1:M
        if (time <= timeend)
            Vxnew(c,N) = Vxnew(c,N-1);
            Vynew(c,N) = Vynew(c,N-1); 
            Cnew(c, N) = Cnew(c, N-1);
            Nnew(c, N) = Nnew(c, N-1);
        end;         
    end;
    
    for i = 2:1:M-1 
        for j = 2:1:N-1  
            if (granitsy(i,j,dx,L,H,img) == 1)
                if (time <= timeend) 
                    %������� �����
                    if (granitsy(i+1,j,dx,L,H,img) == 0)
                        Cnew(i,j) = (-1/Prd*(Cnew(i+1,j)-Cnew(i-1,j))/(2*dx)-Hr*Nnew(i,j))/(-K*(1-Nnew(i,j)/Nx));
                        Vxnew(i,j) = l*(Vxnew(i+1,j)-Vxnew(i-1,j))/(2*dx);
                        %Vynew(i,j) = l*(Vynew(i+1,j)-Vynew(i-1,j))/(2*dx);
                    end;
                    %������� ������
                    if (granitsy(i-1,j,dx,L,H,img) == 0)
                        Cnew(i,j) = (-1)*(-1/Prd*(Cnew(i+1,j)-Cnew(i-1,j))/(2*dx)-Hr*Nnew(i,j))/(-K*(1-Nnew(i,j)/Nx));
                        Vxnew(i,j) = (-1)*l*(Vxnew(i+1,j)-Vxnew(i-1,j))/(2*dx);
                        %Vynew(i,j) = (-1)*l*(Vynew(i+1,j)-Vynew(i-1,j))/(2*dx);
                    end;
                end;
            end;
        end;
    end;
    
    
    
%calculating pressure and velocitycorrection
    for nIter = 1:1:2
        for i = 2:1:M-1
            for j = 2:1:N-1
                if (granitsy(i,j,dx,L,H,img) == 1)
                    if (time <= timeend)
                        DeltaP(i,j) = -gamma*dt/(2*dx)*(Vxold(i+1,j)-Vxold(i-1,j)+Vyold(i,j+1)-Vyold(i,j-1)); 
                    end;                    
                end;
            end;
        end;
        
%calculating correction velocity
        for i=2:1:M-1
            for j=2:1:N-1
                if (granitsy(i,j,dx,L,H,img) == 1)
                    if (time <= timeend)
                        Vxnew(i,j) = Vxnew(i,j)-1/Re*dt*P/(M0*V/R)*(DeltaP(i+1,j)-DeltaP(i-1,j))/(2*dx);
                        Vynew(i,j) = Vynew(i,j)-1/Re*dt*P/(M0*V/R)*(DeltaP(i,j+1)-DeltaP(i,j-1))/(2*dx);
                    end;                    
                end;
            end;
        end;
    end;
   
    
    
    for i=1:1:M
        for j=1:1:N
            if (time <= timeend)
                Vxold(i,j) = Vxnew(i,j);
                Vyold(i,j) = Vynew(i,j);
                Pold(i,j) = Pnew(i,j);
                Cold(i,j) = Cnew(i,j);
                Nold(i,j) = Nnew(i,j);
            end;            
        end;
    end;               
 
      
    time = time + dt;    
    k = k+1;    
    if (time > timeend) flagexit = 0; end;
    
% %     if ((time - timek) > 0.01)
% %         for t = 1:1:N
% %             mas(t) = sqrt(Vyold(M/2,t)^2+Vxold(M/2,t)^2);            
% %         end;
% %         t = 1:1:N;    
% %         h = figure(4);   
% %         clf    
% %         plot(t, mas,'-b');
% %         grid on;
% %         title(strcat('V iteration =  ', num2str(k)));
% %         xlabel('H'); 
% %         ylabel('V');  
% %     end;
    
    if (time <= timeend)
    if ((time - timek) > 0.0001)
        timek = time;       
        h = figure(1);   
        clf    
        subplot(2,2,1);
        surf(X,Y,Vxold);
        title(strcat('Vx time =  ', num2str(time)));
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading interp
        shading flat
        view([0 90]);
        colorbar
    
        subplot(2,2,2)
        surf(X,Y,Vyold);
        shading interp
        title(strcat('Vy iteration =  ', num2str(k)));
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading flat
        view([0 90]);
        colorbar
    
        subplot(2,2,3)
        surf(X,Y,Pold);
        title('P');
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading interp
        shading flat
        view([0 90]);
        colorbar
    
        subplot(2,2,4)
        set(gcf(), 'Renderer', 'painters')
        surf(X,Y,Cold);
        title('C');
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading interp
        shading flat
        view([0 90]);
        colorbar
        
        h = figure(2);   
        clf            
        surf(X,Y,Nold);
        title(strcat('N time =  ', num2str(time)));
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading interp
        shading flat
        view([0 90]);
        colorbar
    end;
    end;
    
end;


time;