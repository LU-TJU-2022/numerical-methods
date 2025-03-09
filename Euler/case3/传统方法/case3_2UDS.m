% clc;
clear;
tic;
X=12800;%//length of the domain
Y=X;
T=12000;%//time of simulation
Nx=128;%//the number of cells in space
Ny=Nx;
M=Nx/128*120000;%//the number of cells in time,time step
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
c_exact=zeros(Nx+1,Ny+1);
c2=zeros(Nx+1,Ny+1);
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
                Cx_i_plu_1_2=3/2*c1(i,j)-1/2*c1(i-1,j);
            else 
                Cx_i_plu_1_2=3/2*c1(i+1,j)-1/2*c1(i+2,j);
            end
            
            if ax_i_min_1_2>0
                Cx_i_min_1_2=3/2*c1(i,j)-1/2*c1(i+1,j);
            else
                Cx_i_min_1_2=3/2*c1(i-1,j)-1/2*c1(i-2,j);
            end
            
            if ay_j_plu_1_2>0
                Cy_j_plu_1_2=3/2*c1(i,j)-1/2*c1(i,j-1);
            else 
                Cy_j_plu_1_2=3/2*c1(i,j+1)-1/2*c1(i,j+2);
            end
            
            if ay_j_min_1_2>0
                Cy_j_min_1_2=3/2*c1(i,j)-1/2*c1(i,j+1);
            else
                Cy_j_min_1_2=3/2*c1(i,j-1)-1/2*c1(i,j-2);
            end
        
        Fx_i_plu_1_2=1/2*(u(i,j)+u(i+1,j))*Cx_i_plu_1_2;
        Fx_i_min_1_2=1/2*(u(i-1,j)+u(i,j))*Cx_i_min_1_2;
        Fy_j_plu_1_2=1/2*(v(i,j)+v(i,j+1))*Cy_j_plu_1_2;
        Fy_j_min_1_2=1/2*(v(i,j-1)+v(i,j))*Cy_j_min_1_2;
   %%%%%%%%%%%%%%% ������ɢ��
        diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
        diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
 %%%%%%%%%%%%%
%       c1(i,j)=c1(i,j)-(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*( Fy_i_plu_1_2-Fy_i_min_1_2)+(dt/(dx^2))*D*diffu_x+(dt/(dy^2))*D*diffu_y;
        c2(i,j)=c1(i,j)-(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*( Fy_j_plu_1_2-Fy_j_min_1_2)+2*(dt/(dx^2))*D*diffu_x+2*(dt/(dy^2))*D*diffu_y;
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


close all;    % �ر�����ͼ�δ���
figure;       % ����һ���µ�ͼ�δ���
contourf(x1,y1,c1');
saveas(gcf,"2UDS1.svg", 'svg');

close all;    % �ر�����ͼ�δ���
figure;       % ����һ���µ�ͼ�δ���
hold on;
box on;
plot(x1,c_exact(:,Ny/4*3+1),':k');
plot(x1,d(:,Ny/4*3+1),'-k');
plot(x1,c1(:,Ny/4*3+1),'--ko');
hold off
saveas(gcf,"2UDS2.svg", 'svg');

error = 0;
for i=1:Nx+1
    for j=1:Ny+1
        error = error + abs(c1(i,j)-d(i,j))^2;
    end
end
error = error/((Nx+1)*(Ny+1));
disp(error);
disp(log(error));
disp(sqrt(dx^2+dy^2));
disp(log(sqrt(dx^2+dy^2)));
disp(max(max(c1)));
disp(min(min(c1)));