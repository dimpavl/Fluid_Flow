clc;
clear;
%Карта нанотрубки
img = imread('Card.bmp');
 
%Характерные параметры
V = 100; % скорость
p = 10^3; % плотность
M0 = 10^(-3); %вязкость
R = 100*10^(-9);%масштаб длины 
D = 10^(-6);%коэффициент диффузии
nu = M0/p; %кинематическая вязкость
Re = R*V/nu;%число Рейнольдса
Prd = R*V/D;%Пранкель диффузионный
P = 10^6; % характерное давление 1 атмосфера
T = R/V; %масштаб по времени
%Параметры для уравнения концентрации жидкости
alpha = 1;

Ce = 0.1;
C0 = 1;



m0 = 1;
u = 0;% коэффициент влияния частиц
u1 = 0.5;

P0 = 0.1;%давление на входе относительное
gamma = 10^(floor(log10(P0)));

%параметры сетки
dx =10^(-2); 
H = 1; L = 0.4;
A = 10^(3);

% Временной отрезок
timeend = 7.1;


N = floor(H/dx);%Продольный размер
M = floor(L/dx);%Поперечный размер
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

Vxold1(M, N) = 0;
Vyold1(M, N) = 0;
Pold1(M, N) = 0;
Vxnew1(M, N) = 0;
Vynew1(M, N) = 0;
Pnew1(M, N) = 0; 
Cnew1(M, N) = 0;
Cold1(M, N) = 0;
m1(M, N) = 0;
DeltaP1(M, N) = 0;
[X,Y] = MeshGrid(1:1:N,1:1:M);

fprintf(1,'Начальные значения\n');
fprintf(1,'-----------------------------------------\n');
fprintf(1,'Для системы уравненя на концентрацию\n');
fprintf(1,'-----------------------------------------\n');
fprintf(1,'alpha\n');
alpha
fprintf(1,'Prd пранкель диффузионный\n');
Prd
fprintf(1,'Ce\n');
Ce
fprintf(1,'Начальная концентрация C0\n');
C0
fprintf(1,'-----------------------------------------\n');
fprintf(1,'Для системы уравнений Навье-Стокса\n');
fprintf(1,'-----------------------------------------\n');
fprintf(1,'Коэффициент кинематической вязкости(начальный)\n');
m0

fprintf(1,'Шаг сетки\n');
dx
fprintf(1,'Давление на входе трубки\n');
P0
fprintf(1,'Продольный размер\n');
H
fprintf(1,'Поперечный размер\n');
L
fprintf(1,'gamma \n');
gamma
fprintf(1,'-----------------------------------------\n');
fprintf(1,'Временной промежуток\n');
[0 timeend]
fprintf(1,'-----------------------------------------\n');
if (u == 0)
    fprintf(1,'Постоянная вязкость\n');
else
    fprintf(1,'Переменная вязкость\n');
end;
%Начальные условия
   for i = 1:1:M     
      if (granitsy(i,1,dx,L,H,img) == 1) 
        Vynew(i,1) = A*((i-1)*dx)^(2)*(L-(i)*dx)^(2);
        Vxnew(i,1)=0;
        Pnew(i,1) = P0;
        Cnew(i,1) = C0;
        Vynew1(i,1) = A*((i-1)*dx)^(2)*(L-(i)*dx)^(2);
        Vxnew1(i,1)=0;
        Pnew1(i,1) = P0;
        Cnew1(i,1) = C0;
      end;
      if (granitsy(i,1,dx,L,H,img) == 1)
        Pnew(i,N) = 0;
        Pnew1(i,N) = 0;
      end;      
   end;
    
  for i = 1:1:M 
        for j = 1:1:N
            m(i,j) = m0;
            m1(i,j) = m0;
        end;
  end;

  
time = 0;
time1 = 0;
timek = 0;
timek1 = 0;
k = 0;
flagexit = 1;
while (flagexit == 1)    
    if (time > timeend && time1 > timeend) flagexit = 0; end;
    if (flagexit == 0) continue; end;
    
    maxM = -Inf; maxM1 = -Inf;
    maxV = -Inf; maxV1 = -Inf;
    maxC = -Inf; maxC1 = -Inf;
    maxP = -Inf; maxP1 = -Inf;
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
               end;
               if (time1 <= timeend)
                    m1(i,j) = m0*(1+u1*Cold1(i,j));                
                    t1 = sqrt(Vxold1(i,j)^2+Vyold1(i,j)^2);                
                    if (t1 > maxV1) maxV1 = t1; end;
                    if (m1(i,j) > maxM1) maxM1 = m1(i,j); end;
                    if (Cold1(i,j) > maxC1) maxC1 = Cold1(i,j); end;
                    if (Pold1(i,j) > maxP1) maxP1 = Pold1(i,j); end;
               end;
           end;
        end;
    end;
    
    if (time <= timeend)
        dt=0.1*min([dx/maxV, dx*dx/(nu)*0.1, dx^2/(1/Prd*maxC)]);
    end;
    if (time1 <= timeend)
        dt1=0.1*min([dx/maxV1, dx*dx/(nu)*0.1, dx^2/(1/Prd*maxC1)]);
    end; 
    
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
                end;
                if (time1 <= timeend)
                    Vxnew1(i,j) =  dt1*(-0.5*(Vxold1(i,j)+abs(Vxold1(i,j))*(Vxold1(i,j)-Vxold1(i-1,j))/dx...
                                -0.5*(Vxold1(i,j)-abs(Vxold1(i,j)))*(Vxold1(i+1,j)-Vxold1(i,j))/dx...
                                -0.5*(Vyold1(i,j)+abs(Vyold1(i,j)))*(Vxold1(i,j)-Vxold1(i,j-1))/dx...
                                -0.5*(Vyold1(i,j)-abs(Vyold1(i,j)))*(Vxold1(i,j+1)-Vxold1(i,j))/dx))...             
                                +1/Re*m1(i,j)*dt1*(Vxold1(i+1,j)+Vxold1(i-1,j)+Vxold1(i,j+1)+Vxold1(i,j-1)-4*Vxold1(i,j))/(dx)^2....
                                -1/Re*dt1*P/(M0*V/R)*(Pold1(i+1,j)-Pold1(i-1,j))/(2*dx)...
                                +1/Re*2*dt1*(m1(i+1,j)-m1(i-1,j))/(2*dx)*(Vxold1(i+1,j)-Vxold1(i-1,j))/(2*dx)...
                                +1/Re*dt1*(m1(i,j+1)-m1(i,j-1))/(2*dx)*(Vyold1(i+1,j)-Vyold1(i-1,j)+Vxold1(i,j+1)-Vxold1(i,j-1))/(2*dx)...    
                                +Vxold1(i,j);
  
                     Vynew1(i,j) = dt1*(-0.5*(Vxold1(i,j)+abs(Vxold1(i,j)))*(Vyold1(i,j)-Vyold1(i-1,j))/dx...
                                -0.5*(Vxold1(i,j)-abs(Vxold1(i,j)))*(Vyold1(i+1,j)-Vyold1(i,j))/dx...
                                -0.5*(Vyold1(i,j)+abs(Vyold1(i,j)))*(Vyold1(i,j)-Vyold1(i,j-1))/dx...
                                -0.5*(Vyold1(i,j)-abs(Vyold1(i,j)))*(Vyold1(i,j+1)-Vyold1(i,j))/dx)...
                                +1/Re*m1(i,j)/p*dt1/(dx)^2*(Vyold1(i+1,j)+Vyold1(i-1,j)+Vyold1(i,j+1)+Vyold1(i,j-1)-4*Vyold1(i,j))...
                                -1/Re*dt1*P/(M0*V/R)*(Pold1(i,j+1)-Pold1(i,j-1))/(2*dx)...
                                +1/Re*dt1*(m1(i+1,j)-m1(i-1,j))/(2*dx)*(Vyold1(i+1,j)-Vyold1(i-1,j)+Vxold1(i,j+1)-Vxold1(i,j-1))/(2*dx)...
                                +1/Re*2*dt1*(m1(i,j+1)-m1(i,j-1))/(2*dx)*(Vyold1(i,j+1)-Vyold1(i,j-1))/(2*dx)...
                                +Vyold1(i,j);
                
                    Pnew1(i,j) = -gamma*dt1/(2*dx)*(Vxold1(i+1,j)-Vxold1(i-1,j)+Vyold1(i,j+1)-Vyold1(i,j-1))+Pold1(i,j);
                
                    Cnew1(i,j) = dt1*1/Prd*(Cold1(i+1,j)+Cold1(i-1,j)+Cold1(i,j+1)+Cold1(i,j-1)-4*Cold1(i,j))/(dx^2)...
                                +alpha*dt1*(Ce-Cold1(i,j))+Cold1(i,j)+dt1*...
                                (-0.5*(Vxold1(i,j)+abs(Vxold1(i,j))*(Cold1(i,j)-Cold1(i-1,j))/dx...
                                -0.5*(Vxold1(i,j)-abs(Vxold1(i,j)))*(Cold1(i+1,j)-Cold1(i,j))/dx...
                                -0.5*(Vyold1(i,j)+abs(Vyold1(i,j)))*(Cold1(i,j)-Cold1(i,j-1))/dx...
                                -0.5*(Vyold1(i,j)-abs(Vyold1(i,j)))*(Cold1(i,j+1)-Cold1(i,j))/dx));
                end;
            end;
        end;
    end;
    for c = 1:1:M
        if (time <= timeend)
            Vxnew(c,N) = Vxnew(c,N-1);
            Vynew(c,N) = Vynew(c,N-1); 
            Cnew(c, N) = Cnew(c, N-1);
        end;
         if (time1 <= timeend)
            Vxnew1(c,N) = Vxnew1(c,N-1);
            Vynew1(c,N) = Vynew1(c,N-1); 
            Cnew1(c, N) = Cnew1(c, N-1);
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
                    if (time1<=timeend)
                        DeltaP1(i,j) = -gamma*dt1/(2*dx)*(Vxold1(i+1,j)-Vxold1(i-1,j)+Vyold1(i,j+1)-Vyold1(i,j-1));
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
                    if (time1 <= timeend)
                        Vxnew1(i,j) = Vxnew1(i,j)-1/Re*dt1*P/(M0*V/R)*(DeltaP1(i+1,j)-DeltaP1(i-1,j))/(2*dx);
                        Vynew1(i,j) = Vynew1(i,j)-1/Re*dt1*P/(M0*V/R)*(DeltaP1(i,j+1)-DeltaP1(i,j-1))/(2*dx);
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
            end;
            if (time1 <= timeend)
                Vxold1(i,j) = Vxnew1(i,j);
                Vyold1(i,j) = Vynew1(i,j);
                Pold1(i,j) = Pnew1(i,j);
                Cold1(i,j) = Cnew1(i,j);
            end;
        end;
    end;               
 
      
    time = time + dt;
    time1 = time1+dt1;
    k = k+1;    
    if (time > timeend && time1 > timeend) flagexit = 0; end;
    
    if ((time - timek) > 0.01 || (time1 - timek1) > 0.01)
        for t = 1:1:N
            mas(t) = sqrt(Vyold(M/2,t)^2+Vxold(M/2,t)^2);
            mas1(t) = sqrt(Vyold1(M/2,t)^2+Vxold1(M/2,t)^2);
        end;
        t = 1:1:N;    
        h = figure(4);   
        clf    
        plot(t, mas,'-b', t, mas1,'-r');
        grid on;
        title(strcat('V iteration =  ', num2str(k)));
        xlabel('H'); 
        ylabel('V');  
    end;
    
    if (time <= timeend)
    if ((time - timek) > 0.01)
        timek = time;       
        h = figure(5);   
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
    end;
    end;
    
    if (time1 <= timeend)
    if ((time1 - timek1) > 0.01)
        timek1 = time1;       
        h = figure(6);   
        clf    
        subplot(2,2,1);
        surf(X,Y,Vxold1);
        title(strcat('Vx time =  ', num2str(time1)));
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading interp
        shading flat
        view([0 90]);
        colorbar
    
        subplot(2,2,2)
        surf(X,Y,Vyold1);
        shading interp
        title(strcat('Vy iteration =  ', num2str(k)));
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading flat
        view([0 90]);
        colorbar
    
        subplot(2,2,3)
        surf(X,Y,Pold1);
        title('P');
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading interp
        shading flat
        view([0 90]);
        colorbar
    
        subplot(2,2,4)
        set(gcf(), 'Renderer', 'painters')
        surf(X,Y,Cold1);
        title('C');
        xlabel( strcat('H= ', num2str(H)) ); ylabel( strcat('L= ', num2str(L)) );
        shading interp
        shading flat
        view([0 90]);
        colorbar       
    end;
    end;
    
  
end;


time;

