%------WENO:Weighted Essentially Non-Oscillatory--------%
%-----Jiang & Shu 1996----%
%-----5阶WENO-------%
clear
clc


tic;

d1=1/10;
d2=3/5;
d3=3/10;
e=1e-6;


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

x0=-1.5;
y0=-1.5;
% t=dt*(1:M+1);

a=1.5;

x1=zeros(1,Nx+1);
y1=zeros(Ny+1,1);
c1=zeros(Nx+1,Ny+1);
c2=zeros(Nx+1,Ny+1);
u=zeros(Nx+1,Ny+1)+0.1;
v=zeros(Nx+1,Ny+1)+0.1;

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

% figure(1)
% contourf(x1,y1,c1')
% mesh(x1,y1,c1')

for n=1:M+1
    disp(n);
    for i=4:Nx-1
        for j=4:Ny-1
            %------Fxi+1/2------%
            beta_x1=13/12*(u(i-2,j)*c1(i-2,j)-2*u(i-1,j)*c1(i-1,j)+u(i,j)*c1(i,j)).^2+1/4*(u(i-2,j)*c1(i-2,j)-4*u(i-1,j)*c1(i-1,j)+3*u(i,j)*c1(i,j)).^2;
            beta_x2=13/12*(u(i-1,j)*c1(i-1,j)-2*u(i,j)*c1(i,j)+u(i+1,j)*c1(i+1,j)).^2+1/4*(u(i-1,j)*c1(i-1,j)-u(i+1,j)*c1(i+1,j)).^2;
            beta_x3=13/12*(u(i,j)*c1(i,j)-2*u(i+1,j)*c1(i+1,j)+u(i+2,j)*c1(i+2,j)).^2+1/4*(3*u(i,j)*c1(i,j)-4*u(i+1,j)*c1(i+1,j)+u(i+2,j)*c1(i+2,j)).^2;
            ax1=d1/(e+beta_x1)^2;
            ax2=d2/(e+beta_x2)^2;
            ax3=d3/(e+beta_x3)^2;
            wx1=ax1/(ax1+ax2+ax3);
            wx2=ax2/(ax1+ax2+ax3);
            wx3=ax3/(ax1+ax2+ax3);
            
            Fx_i_plu_1_2_1=1/3*u(i-2,j)*c1(i-2,j)-7/6*u(i-1,j)*c1(i-1,j)+11/6*u(i,j)*c1(i,j);
            Fx_i_plu_1_2_2=-1/6*u(i-1,j)*c1(i-1,j)+5/6*u(i,j)*c1(i,j)+1/3*u(i+1,j)*c1(i+1,j);
            Fx_i_plu_1_2_3=1/3*u(i,j)*c1(i,j)+5/6*u(i+1,j)*c1(i+1,j)-1/6*u(i+2,j)*c1(i+2,j);

            Fx_i_plu_1_2=wx1*Fx_i_plu_1_2_1+wx2*Fx_i_plu_1_2_2+wx3*Fx_i_plu_1_2_3;
            
             %-------Fxi-1/2---------%
            beta_x4=13/12*(u(i-3,j)*c1(i-3,j)-2*u(i-2,j)*c1(i-2,j)+u(i-1,j)*c1(i-1,j)).^2+1/4*(u(i-3,j)*c1(i-3,j)-4*u(i-2,j)*c1(i-2,j)+3*u(i-1,j)*c1(i-1,j)).^2;
            beta_x5=13/12*(u(i-2,j)*c1(i-2,j)-2*u(i-1,j)*c1(i-1,j)+u(i,j)*c1(i,j)).^2+1/4*(u(i-2,j)*c1(i-2,j)-u(i,j)*c1(i,j)).^2;
            beta_x6=13/12*(u(i-1,j)*c1(i-1,j)-2*u(i,j)*c1(i,j)+u(i+1,j)*c1(i+1,j)).^2+1/4*(3*u(i-1,j)*c1(i-1,j)-4*u(i,j)*c1(i,j)+u(i+1,j)*c1(i+1,j)).^2;
            ax4=d1/(e+beta_x4)^2;
            ax5=d2/(e+beta_x5)^2;
            ax6=d3/(e+beta_x6)^2;
            wx4=ax4/(ax4+ax5+ax6);
            wx5=ax5/(ax4+ax5+ax6);
            wx6=ax6/(ax4+ax5+ax6);
            
            Fx_i_min_1_2_1=1/3*u(i-3,j)*c1(i-3,j)-7/6*u(i-2,j)*c1(i-2,j)+11/6*u(i-1,j)*c1(i-1,j);
            Fx_i_min_1_2_2=-1/6*u(i-2,j)*c1(i-2,j)+5/6*u(i-1,j)*c1(i-1,j)+1/3*u(i,j)*c1(i,j);
            Fx_i_min_1_2_3=1/3*u(i-1,j)*c1(i-1,j)+5/6*u(i,j)*c1(i,j)-1/6*u(i+1,j)*c1(i+1,j);
            
            Fx_i_min_1_2=wx4*Fx_i_min_1_2_1+wx5*Fx_i_min_1_2_2+wx6*Fx_i_min_1_2_3;
           
             %--------Fyi+1/2-------%
            beta_y1=13/12*(v(i,j-2)*c1(i,j-2)-2*v(i,j-1)*c1(i,j-1)+v(i,j)*c1(i,j)).^2+1/4*(v(i,j-2)*c1(i,j-2)-4*v(i,j-1)*c1(i,j-1)+3*v(i,j)*c1(i,j)).^2;
            beta_y2=13/12*(v(i,j-1)*c1(i,j-1)-2*v(i,j)*c1(i,j)+v(i,j+1)*c1(i,j+1)).^2+1/4*(v(i,j-1)*c1(i,j-1)-v(i,j+1)*c1(i,j+1)).^2;
            beta_y3=13/12*(v(i,j)*c1(i,j)-2*v(i,j+1)*c1(i,j+1)+v(i,j+2)*c1(i,j+2)).^2+1/4*(3*v(i,j)*c1(i,j)-4*v(i,j+1)*c1(i,j+1)+v(i,j+2)*c1(i,j+2)).^2;
            ay1=d1/(e+beta_y1)^2;
            ay2=d2/(e+beta_y2)^2;
            ay3=d3/(e+beta_y3)^2;
            wy1=ay1/(ay1+ay2+ay3);
            wy2=ay2/(ay1+ay2+ay3);
            wy3=ay3/(ay1+ay2+ay3);
            
            Fy_i_plu_1_2_1=1/3*v(i,j-2)*c1(i,j-2)-7/6*v(i,j-1)*c1(i,j-1)+11/6*v(i,j)*c1(i,j);
            Fy_i_plu_1_2_2=-1/6*v(i,j-1)*c1(i,j-1)+5/6*v(i,j)*c1(i,j)+1/3*v(i,j+1)*c1(i,j+1);
            Fy_i_plu_1_2_3=1/3*v(i,j)*c1(i,j)+5/6*v(i,j+1)*c1(i,j+1)-1/6*v(i,j+2)*c1(i,j+2);
            
            Fy_i_plu_1_2=wy1*Fy_i_plu_1_2_1+wy2*Fy_i_plu_1_2_2+wy3*Fy_i_plu_1_2_3;
            
           
            %--------Fyi-1/2---------%
            beta_y4=13/12*(v(i,j-3)*c1(i,j-3)-2*v(i,j-2)*c1(i,j-2)+v(i,j-1)*c1(i,j-1)).^2+1/4*(v(i,j-3)*c1(i,j-3)-4*v(i,j-2)*c1(i,j-2)+3*v(i,j-1)*c1(i,j-1)).^2;
            beta_y5=13/12*(v(i,j-2)*c1(i,j-2)-2*v(i,j-1)*c1(i,j-1)+v(i,j)*c1(i,j)).^2+1/4*(v(i,j-2)*c1(i,j-2)-v(i,j)*c1(i,j)).^2;
            beta_y6=13/12*(v(i,j-1)*c1(i,j-1)-2*v(i,j)*c1(i,j)+v(i,j+1)*c1(i,j+1)).^2+1/4*(3*v(i,j-1)*c1(i,j-1)-4*v(i,j)*c1(i,j)+v(i,j+1)*c1(i,j+1)).^2;
            ay4=d1/(e+beta_y4)^2;
            ay5=d2/(e+beta_y5)^2;
            ay6=d3/(e+beta_y6)^2;
            wy4=ay4/(ay4+ay5+ay6);
            wy5=ay5/(ay4+ay5+ay6);
            wy6=ay6/(ay4+ay5+ay6);
            
            Fy_i_min_1_2_1=1/3*v(i,j-3)*c1(i,j-3)-7/6*v(i,j-2)*c1(i,j-2)+11/6*v(i,j-1)*c1(i,j-1);
            Fy_i_min_1_2_2=-1/6*v(i,j-2)*c1(i,j-2)+5/6*v(i,j-1)*c1(i,j-1)+1/3*v(i,j)*c1(i,j);
            Fy_i_min_1_2_3=1/3*v(i,j-1)*c1(i,j-1)+5/6*v(i,j)*c1(i,j)-1/6*v(i,j+1)*c1(i,j+1);
            
            Fy_i_min_1_2=wy4*Fy_i_min_1_2_1+wy5*Fy_i_min_1_2_2+wy6*Fy_i_min_1_2_3;
            
            
            diff_x=c1(i-1,j)-2*c1(i,j)+c1(i+1,j);
            diff_y=c1(i,j-1)-2*c1(i,j)+c1(i,j+1);
            c2(i,j)=c1(i,j)-dt/dx*(Fx_i_plu_1_2-Fx_i_min_1_2)-dt/dy*(Fy_i_plu_1_2-Fy_i_min_1_2)+D*dt/(dx^2)*diff_x+D*dt/(dx^2)*diff_y;
        end
    end
    
    c2(3,:)=c2(4,:);
    c2(2,:)=c2(3,:);
    c2(1,:)=c2(2,:);
    c2(Nx,:)=c2(Nx-1,:);
    c2(Nx+1,:)=c2(Nx,:);
    
    c2(:,3)=c2(:,4);
    c2(:,2)=c2(:,3);
    c2(:,1)=c2(:,2);
    c2(:,Ny)=c2(:,Ny-1);
    c2(:,Ny+1)=c2(:,Ny);
    
    c1=c2;
    
end
toc;
%   figure(1)
%   hold on;
% plot(c1(:,17),'r'); 
figure(6);
contourf(x1,y1,c1')

disp(max(max(c1)));
disp(min(min(c1)));
saveas(gcf, 'WENO.svg', 'svg');
