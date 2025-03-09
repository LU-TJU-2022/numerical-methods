clc;
clear;
tic;
X=6400;%//length of the domain
Y=X;
T=3500; 
Nx=128;
Ny=Nx;
% M=Nx/128*6400;
M=3500;
% w=pi/6000;
dx=X/Nx;% //this is the size of the space step
dy=Y/Ny;
sx=200*1;
sy=200*1;
dt=T/M; %//this is the size of the time step
D=0.1;%%%扩散系数
x1=zeros(1,Nx+1); % array allocation
y1=zeros(1,Ny+1);
c1=zeros(Nx+1,Ny+1);
c_exact=zeros(Nx+1,Ny+1);
d=zeros(Nx+1,Ny+1);
u=zeros(Nx+1,Ny+1);
v=zeros(Nx+1,Ny+1);
r=zeros(Nx+1,Ny+1);
x0=1300;
y0=1300;
% rc=0.1;
for i=1:Nx+1
    for j=1:Ny+1
    x1(i)=dx*(i-1); %//there is no cell 0 in scilab that is why we begin with 1
    y1(j)=dy*(j-1);
    c1(i,j)=exp(-(x1(i)-x0)^2/(2*sx^2)-(y1(j)-y0)^2/(2*sy^2));
    u(i,j)=1;
    v(i,j)=1;
    end
end
c2=c1;
c_exact = c1;
% for i=1:Nx+1
%     for j=1:Ny+1
%     c2(i,j)=exp(-(x1(i)-(x0+u(i,j)*dt))^2/(2*sx^2)-(y1(j)-(y0+v(i,j)*dt))^2/(2*sy^2));
%     end
% end
for i=1:Nx+1
    for j=1:Ny+1
%     d(i,j)=exp(-(x1(i)-(x0+u(i,j)*T))^2/(2*sx^2)-(y1(j)-(y0+v(i,j)*T))^2/(2*sy^2));
    d(i,j)=(sx/sqrt(sx^2+2*D*T))*(sy/sqrt(sy^2+2*D*T))*exp(-(x1(i)-(x0+u(i,j)*T))^2/(2*(sx^2+2*D*T))-(y1(j)-(y0+v(i,j)*T))^2/(2*(sy^2+2*D*T)));
    end
end
% [X,Y]=meshgrid(1:Nx+1,1:Ny+1);%%%%meshgrid,注意矩阵行数对应y，列数对应X，比如3*4列，则y=1:3，x=1:4
% figure(1);
% contourf(x1,y1,c1');
% xlabel("x");
% ylabel("y");
% title("Initial distribution");
% 
% figure(4);
% contourf(x1,y1,d');
% xlabel("x");
% ylabel("y");
% title("Exact solution");

for n=1:M
    for i=3:Nx-1
        for j=3:Ny-1
%%%%%%%%%%%%%%% 计算对流项
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
            
            
            if ax_i_plu_1_2>0
                Cx_i_plu_1_2=c1(i,j)+1/4*(c1(i+1,j)-c1(i-1,j))+1/8*(c1(i+1,j)-2*c1(i,j)+c1(i-1,j));
            else 
                Cx_i_plu_1_2=c1(i+1,j)-1/4*(c1(i+2,j)-c1(i,j))+1/8*(c1(i+2,j)-2*c1(i+1,j)+c1(i,j));
            end
            
            if ax_i_min_1_2>0
                Cx_i_min_1_2=c1(i,j)+1/4*(c1(i-1,j)-c1(i+1,j))+1/8*(c1(i-1,j)-2*c1(i,j)+c1(i+1,j));
            else
                Cx_i_min_1_2=c1(i-1,j)-1/4*(c1(i-2,j)-c1(i,j))+1/8*(c1(i-2,j)-2*c1(i-1,j)+c1(i,j));
            end
            
            if ay_j_plu_1_2>0
                Cy_j_plu_1_2=c1(i,j)+1/4*(c1(i,j+1)-c1(i,j-1))+1/8*(c1(i,j+1)-2*c1(i,j)+c1(i,j-1));
            else 
                Cy_j_plu_1_2=c1(i,j+1)-1/4*(c1(i,j+2)-c1(i,j))+1/8*(c1(i,j+2)-2*c1(i,j+1)+c1(i,j));
            end
            
            if ay_j_min_1_2>0
                Cy_j_min_1_2=c1(i,j)+1/4*(c1(i,j-1)-c1(i,j+1))+1/8*(c1(i,j-1)-2*c1(i,j)+c1(i,j+1));
            else
                Cy_j_min_1_2=c1(i,j-1)-1/4*(c1(i,j-2)-c1(i,j))+1/8*(c1(i,j-2)-2*c1(i,j-1)+c1(i,j));
            end
%         Cx_i_plu_1_2=1/2*(c2(i,j)+c2(i+1,j));
%         Cx_i_min_1_2=1/2*(c2(i-1,j)+c2(i,j));
%         Cy_j_plu_1_2=1/2*(c2(i,j)+c2(i,j+1));
%         Cy_j_min_1_2=1/2*(c2(i,j-1)+c2(i,j));
        
        Fx_i_plu_1_2=1/2*(u(i,j)+u(i+1,j))*Cx_i_plu_1_2;
        Fx_i_min_1_2=1/2*(u(i-1,j)+u(i,j))*Cx_i_min_1_2;
        Fy_j_plu_1_2=1/2*(v(i,j)+v(i,j+1))*Cy_j_plu_1_2;
        Fy_j_min_1_2=1/2*(v(i,j-1)+v(i,j))*Cy_j_min_1_2;
   %%%%%%%%%%%%%%% 计算扩散项
        diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
        diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
 %%%%%%%%%%%%%
%       c1(i,j)=c1(i,j)-(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*( Fy_i_plu_1_2-Fy_i_min_1_2)+(dt/(dx^2))*D*diffu_x+(dt/(dy^2))*D*diffu_y;
        c2(i,j)=c1(i,j)-(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*( Fy_j_plu_1_2-Fy_j_min_1_2)+2*(dt/(dx^2))*D*diffu_x+2*(dt/(dy^2))*D*diffu_y;
       end
    end
    c1=c2;
%     figure(10);
%     xlabel("x");
%     ylabel("y");
%     hold on
%     if mod(n,200)==0
%         contourf(x1,y1,c3');
%         title(n+"/"+M);
%     end
%     hold off
    
end

 toc;
close all;    % 关闭所有图形窗口
figure;       % 创建一个新的图形窗口
contourf(x1,y1,c1');
saveas(gcf,"QUICK1.svg", 'svg');

close all;    % 关闭所有图形窗口
figure;       % 创建一个新的图形窗口
hold on;
box on;
plot(x1,c_exact(:,27),':k');
plot(x1,d(:,Ny/4*3+1),'-k');
plot(x1,c1(:,Ny/4*3+1),'--ko');
hold off
saveas(gcf,"QUICK2.svg", 'svg');
disp(max(max(c1)));
disp(min(min(c1)));