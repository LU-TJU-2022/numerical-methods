% clc;
clear;
tic;
X=12800;%//length of the domain
Y=X;
T=12000;%//time of simulation
Nx=128;%//the number of cells in space
Ny=Nx;
M=Nx/128*12000;%//the number of cells in time,time step
w=pi/6000;
dx=X/Nx;% //this is the size of the space step
dy=Y/Ny;
sx=141.43*5;
sy=100*5;
dt=T/M; %//this is the size of the time step
D=1;
x1=zeros(1,Nx+1); % array allocation
y1=zeros(1,Ny+1);
c1=zeros(Nx+1,Ny+1);
c2=zeros(Nx+1,Ny+1);
c_exact=zeros(Nx+1,Ny+1);
d=zeros(Nx+1,Ny+1);
u=zeros(Nx+1,Ny+1);
v=zeros(Nx+1,Ny+1);
r=zeros(Nx+1,Ny+1);
x0=0;
y0=3200;

for i=1:Nx+1
    for j=1:Ny+1
        x1(i)=dx*(i-(Nx/2+1));
        y1(j)=dy*(j-(Ny/2+1));
        c1(i,j)=exp(-(x1(i)-x0)^2/(2*sx^2)-(y1(j)-y0)^2/(2*sy^2));
        u(i,j)=-w*y1(j);
        v(i,j)=w*x1(i);
    end
end
c_exact = c1;

for i=1:Nx+1
    for j=1:Ny+1
    d(i,j)=(sx/sqrt(sx^2+2*D*T))*(sy/sqrt(sy^2+2*D*T))*exp(-(x1(i)-x0)^2/(2*(sx^2+2*D*T))-(y1(j)-y0)^2/(2*(sy^2+2*D*T)));
    end
end


for n=1:M
    for i=3:Nx-1
        for j=3:Ny-1
             Qx_i_plu_1_2=(dx/dt)*(c1(i+1,j)-c1(i,j));
            Qx_i_min_1_2=(dx/dt)*(c1(i,j)-c1(i-1,j));
            Qy_j_plu_1_2=(dy/dt)*(c1(i,j+1)-c1(i,j));
            Qy_j_min_1_2=(dy/dt)*(c1(i,j)-c1(i,j-1));

            fx_i_plu_1=u(i+1,j)*c1(i+1,j);
            fx_i_min_1=u(i-1,j)*c1(i-1,j);
            fy_j_plu_1=v(i,j+1)*c1(i,j+1);
            fy_j_min_1=v(i,j-1)*c1(i,j-1);
            fij=u(i,j)*c1(i,j);

            Fx_i_plu_1_2=1/2*(fx_i_plu_1+fij-Qx_i_plu_1_2);
            Fx_i_min_1_2=1/2*(fij+fx_i_min_1-Qx_i_min_1_2);
            Fy_j_plu_1_2=1/2*(fy_j_plu_1+fij-Qy_j_plu_1_2);
            Fy_j_min_1_2=1/2*(fij+fy_j_min_1-Qy_j_min_1_2);
            
            diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
            diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
 %%%%%%%%%%%%%
%         c2(i,j)=c1(i,j)-1/2*(c1(i+1,j)+c1(i-1,j))-1/2*(c1(i,j+1)+c1(i,j-1))+dt/(2*dx)*(u(i+1,j)*c1(i+1,j)-u(i-1,j)*c1(i-1,j))+dt/(2*dy)*(v(i,j+1)*c1(i,j+1)-v(i,j-1)*c1(i,j-1));
        c2(i,j)=-c1(i,j)+1/2*(c1(i+1,j)+c1(i-1,j))+1/2*(c1(i,j+1)+c1(i,j-1))-(u(i,j)*dt)/(2*dx)*(c1(i+1,j)-c1(i-1,j))-(v(i,j)*dt)/(2*dy)*(c1(i,j+1)-c1(i,j-1))+(dt/(dx^2))*D*diffu_x+(dt/(dy^2))*D*diffu_y;
        end
    end
    c1=c2;
%     if mod(n,100)==0
%         figure(11)
%         contourf(x1,y1,c1');
%     end
end
%%%%%%%%%%  Pass the C value of point d to node(i,j)
toc;

close all;    % 关闭所有图形窗口
figure;       % 创建一个新的图形窗口
contourf(x1,y1,c1');
saveas(gcf,"LFS1.svg", 'svg');

close all;    % 关闭所有图形窗口
figure;       % 创建一个新的图形窗口
hold on;
box on;
plot(x1,c_exact(:,Ny/4*3+1),':k');
plot(x1,d(:,Ny/4*3+1),'-k');
plot(x1,c1(:,Ny/4*3+1),'--ko');
hold off
saveas(gcf,"LFS2.svg", 'svg');
error = 0;
for i=1:Nx+1
    for j=1:Ny+1
        error = error + abs(c1(i,j)-d(i,j))^2;
    end
end
error = error/((Nx+1)*(Ny+1));

disp(max(max(c1)));
disp(min(min(c1)));
disp(error);