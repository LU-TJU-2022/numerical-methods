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

for n=1:M
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
             if c1(i+1,j)~=c1(i,j)
                ax_i_plu_1_2=(u(i+1,j)*c1(i+1,j)-u(i,j)*c1(i,j))/(c1(i+1,j)-c1(i,j));
            else
                ax_i_plu_1_2=u(i,j)*c1(i,j);
            end
            if c1(i-1,j)~=c1(i,j)
                ax_i_min_1_2=(u(i,j)*c1(i,j)-u(i-1,j)*c1(i-1,j))/(c1(i,j)-c1(i-1,j));
            else
                ax_i_min_1_2=u(i,j)*c1(i,j); 
            end
            if c1(i,j+1)~=c1(i,j)
                ay_j_plu_1_2=(v(i,j+1)*c1(i,j+1)-v(i,j)*c1(i,j))/(c1(i,j+1)-c1(i,j));
            else
                ay_j_plu_1_2=v(i,j)*c1(i,j);
            end
            if c1(i,j)~=c1(i,j-1)
                ay_j_min_1_2=(v(i,j)*c1(i,j)-v(i,j-1)*c1(i,j-1))/(c1(i,j)-c1(i,j-1));
            else
                ay_j_min_1_2=v(i,j)*c1(i,j);
            end
   %%%%%%%%%%%%%%% 计算对流项
        Qx_i_plu_1_2=abs(ax_i_plu_1_2)*(c1(i+1,j)-c1(i,j));
        Qx_i_min_1_2=abs(ax_i_min_1_2)*(c1(i,j)-c1(i-1,j));
        Qy_j_plu_1_2=abs(ay_j_plu_1_2)*(c1(i,j+1)-c1(i,j));
        Qy_j_min_1_2=abs(ay_j_min_1_2)*(c1(i,j)-c1(i,j-1));
        
        fx_i_plu_1=u(i+1,j)*c1(i+1,j);
        fx_i_min_1=u(i-1,j)*c1(i-1,j);
        fy_j_plu_1=v(i,j+1)*c1(i,j+1);
        fy_j_min_1=v(i,j-1)*c1(i,j-1);
        fij=u(i,j)*c1(i,j);
        
        Fx_i_plu_1_2=1/2*(fx_i_plu_1+fij-Qx_i_plu_1_2);
        Fx_i_min_1_2=1/2*(fij+fx_i_min_1-Qx_i_min_1_2);
        Fy_j_plu_1_2=1/2*(fy_j_plu_1+fij-Qy_j_plu_1_2);
        Fy_j_min_1_2=1/2*(fij+fy_j_min_1-Qy_j_min_1_2);

   %%%%%%%%%%%%%%% 计算扩散项
        diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
        diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
 %%%%%%%%%%%%%
%       c2(i,j)=c1(i,j)-(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*( Fy_i_plu_1_2-Fy_i_min_1_2)+(dt/(dx^2))*D*diffu_x+(dt/(dy^2))*D*diffu_y;
%         c2(i,j)=c1(i,j)-2*(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-2*(dt/dy)*( Fy_j_plu_1_2-Fy_j_min_1_2)+2*(dt/(dx^2))*D*diffu_x+2*(dt/(dy^2))*D*diffu_y;
%         c2(i,j)=c1(i,j)-dt/(2*dx)*(u(i+1,j)*c1(i+1,j)-u(i-1,j)*c1(i-1,j)-abs(ax_i_plu_1_2)*(c1(i+1,j)-c1(i,j))+abs(ax_i_min_1_2)*(c1(i,j)-c1(i-1,j)))-dt/(2*dy)*(v(i,j+1)*c1(i,j+1)-v(i,j-1)*c1(i,j-1)-abs(ay_j_plu_1_2)*(c1(i,j+1)-c1(i,j))+abs(ay_j_min_1_2)*(c1(i,j)-c1(i,j-1)));
        c2(i,j)=c1(i,j)-dt/(dx)*(Fx_i_plu_1_2-Fx_i_min_1_2)-dt/(dy)*(Fy_j_plu_1_2-Fy_j_min_1_2);
       end
    end
    c1=c2;
end
 
 toc;
contourf(x1,y1,c1')
% 设置坐标轴范围
xlim([-3 3]);  % 设置x轴的范围为-3到3
ylim([-3 3]);  % 设置y轴的范围为-3到3

disp(max(max(c1)));
disp(min(min(c1)));
saveas(gcf, 'UDS.svg', 'svg');

        