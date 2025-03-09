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
D=0;%%%��ɢϵ��
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
[X,Y]=meshgrid(1:Nx+1,1:Ny+1);%%%%meshgrid,ע�����������Ӧy��������ӦX������3*4�У���y=1:3��x=1:4
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
 %%%%%%%%%%%%%%% ���������
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
   %%%%%%%%%%%%%%% ������ɢ��
        diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
        diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
 %%%%%%%%%%%%%
%       c1(i,j)=c1(i,j)-(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*( Fy_i_plu_1_2-Fy_i_min_1_2)+(dt/(dx^2))*D*diffu_x+(dt/(dy^2))*D*diffu_y;
        c2(i,j)=c1(i,j)-(dt/dx)*( Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*( Fy_j_plu_1_2-Fy_j_min_1_2)+2*(dt/(dx^2))*D*diffu_x+2*(dt/(dy^2))*D*diffu_y;
       end
    end
    c1=c2;
end
 
 toc;
contourf(x1,y1,c1')
% ���������᷶Χ
xlim([-3 3]);  % ����x��ķ�ΧΪ-3��3
ylim([-3 3]);  % ����y��ķ�ΧΪ-3��3

disp(max(max(c1)));
disp(min(min(c1)));
saveas(gcf, 'QUICK.svg', 'svg');
        