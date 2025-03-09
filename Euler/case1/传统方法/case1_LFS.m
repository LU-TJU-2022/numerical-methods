clc;
clear;
tic;
X=6;%//length of the domain
Y=6;
T=2.8;
Nx=64;
Ny=64;
M=560;
dx=X/Nx;% //this is the size of the space step
dy=Y/Ny;
dt=T/M; %//this is the size of the time step
D=0;%%%扩散系数
x1=zeros(1,Nx+1); % array allocation
y1=zeros(1,Ny+1);
c1=zeros(Nx+1,Ny+1);
u=zeros(Nx+1,Ny+1);
v=zeros(Nx+1,Ny+1);

x0=-1.5;
y0=-1.5;
% t=dt*(1:M+1);

a=1.5;
r=zeros(Nx+1,Ny+1);
for i=1:Nx+1
    for j=1:Ny+1
    x1(i)=dx*(i-(Nx/2+1)); %//there is no cell 0 in scilab that is why we begin with 1
    y1(j)=dy*(j-(Ny/2+1));
    u(i,j)=1;
    v(i,j)=1;
        if abs(x1(i)-x0)<=a/2 && abs(y1(j)-y0)<=a/2
            c1(i,j)=10;
        else
            c1(i,j)=0;
        end
    end
end

c2=c1;
[X,Y]=meshgrid(1:Nx+1,1:Ny+1);%%%%meshgrid,注意矩阵行数对应y，列数对应X，比如3*4列，则y=1:3，x=1:4
figure(2);
contourf(x1,y1,c1')

for n=1:1
    c1(1,:)=0;
    c1(Nx+1,:)=0;
    c1(:,Ny+1)=0;
    c1(:,1)=0;
    
    c1(2,:)=0;
    c1(Nx,:)=0;
    c1(:,Ny)=0;
    c1(:,2)=0;
    for i=3:Nx-1
        for j=3:Ny-1
%             Qx_i_plu_1_2=(dx/dt)*(c1(i+1,j)-c1(i,j));
%             Qx_i_min_1_2=(dx/dt)*(c1(i,j)-c1(i-1,j));
%             Qy_j_plu_1_2=(dy/dt)*(c1(i,j+1)-c1(i,j));
%             Qy_j_min_1_2=(dy/dt)*(c1(i,j)-c1(i,j-1));
% 
%             fx_i_plu_1=u(i+1,j)*c1(i+1,j);
%             fx_i_min_1=u(i-1,j)*c1(i-1,j);
%             fy_j_plu_1=v(i,j+1)*c1(i,j+1);
%             fy_j_min_1=v(i,j-1)*c1(i,j-1);
%             fij=u(i,j)*c1(i,j);
% 
%             Fx_i_plu_1_2=1/2*(fx_i_plu_1+fij-Qx_i_plu_1_2);
%             Fx_i_min_1_2=1/2*(fij+fx_i_min_1-Qx_i_min_1_2);
%             Fy_j_plu_1_2=1/2*(fy_j_plu_1+fij-Qy_j_plu_1_2);
%             Fy_j_min_1_2=1/2*(fij+fy_j_min_1-Qy_j_min_1_2);
            
            diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
            diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
 %%%%%%%%%%%%%
%         c2(i,j)=c1(i,j)-1/2*(c1(i+1,j)+c1(i-1,j))-1/2*(c1(i,j+1)+c1(i,j-1))+dt/(2*dx)*(u(i+1,j)*c1(i+1,j)-u(i-1,j)*c1(i-1,j))+dt/(2*dy)*(v(i,j+1)*c1(i,j+1)-v(i,j-1)*c1(i,j-1));
        c2(i,j)=-c1(i,j)+1/2*(c1(i+1,j)+c1(i-1,j))+1/2*(c1(i,j+1)+c1(i,j-1))-(u(i,j)*dt)/(2*dx)*(c1(i+1,j)-c1(i-1,j))-(v(i,j)*dt)/(2*dy)*(c1(i,j+1)-c1(i,j-1))+(dt/(dx^2))*D*diffu_x+(dt/(dy^2))*D*diffu_y;
       end
    end
    c1=c2;
    figure(100)
    contourf(x1,y1,c1);
    title(n);
end
 
 toc;
%   figure(1)
%   hold on;
% plot(c1(:,17),'r'); 
figure(6);
contourf(x1,y1,c1')

disp(max(max(c1)));
disp(min(min(c1)));

        