function [Limiter1x, Limiter1y] = TCDF(c1, dx, dy, i, j)
    % TCDF ���������㺯��
    % ����:
    %   c1  - ���������ϵ�Ũ�ȳ�����
    %   dx, dy - ������
    %   i, j - ��ǰ����������
    % ���:
    %   Limiter1x - x ����� TCDF ������ֵ
    %   Limiter1y - y ����� TCDF ������ֵ
    
    XHx_i=(c1(i,j)-c1(i-1,j))/(c1(i+1,j)-c1(i,j));%%%seita(i)
    XHy_i=(c1(i,j)-c1(i,j-1))/(c1(i,j+1)-c1(i,j));%%%seita(i)
    %%%%%%%����limiter����i-1
    XHx_i_min_1=(c1(i-1,j)-c1(i-2,j))/(c1(i,j)-c1(i-1,j));%%%seita(i-1)
    XHy_i_min_1=(c1(i,j-1)-c1(i,j-2))/(c1(i,j)-c1(i,j-1));%%%seita(i-1)
    %%%%%%%%%����limiter����i+1
    XHx_i_plu_1=(c1(i+1,j)-c1(i,j))/(c1(i+2,j)-c1(i+1,j));%%%seita(i+1)
    XHy_i_plu_1=(c1(i,j+1)-c1(i,j))/(c1(i,j+2)-c1(i,j+1));%%%seita(i+1)
    
      if isnan(XHx_i)  %%%  XHx = 0/0
        Limiter1x_i=2*(c1(i+1,j)-c1(i,j))/dx;
      elseif XHx_i==inf
        Limiter1x_i=2*(c1(i+1,j)-c1(i,j))/dx;
      elseif XHx_i==-inf
          Limiter1x_i=1*(c1(i+1,j)-c1(i,j))/dx;
      elseif XHx_i<0
          Limiter1x_i=(XHx_i*XHx_i+XHx_i)/(1+XHx_i*XHx_i)*(c1(i+1,j)-c1(i,j))/dx;
      elseif 0<=XHx_i && XHx_i<=1/2
          Limiter1x_i=(XHx_i*XHx_i*XHx_i-2*XHx_i*XHx_i+2*XHx_i)*(c1(i+1,j)-c1(i,j))/dx;
      elseif 1/2<XHx_i && XHx_i<=2
          Limiter1x_i=(3/4*XHx_i+1/4)*(c1(i+1,j)-c1(i,j))/dx;
      elseif XHx_i>2
          Limiter1x_i=(2*XHx_i*XHx_i-2*XHx_i-9/4)/(XHx_i*XHx_i-XHx_i-1)*(c1(i+1,j)-c1(i,j))/dx;
      end


     if isnan(XHy_i)
         Limiter1y_i=2*(c1(i,j+1)-c1(i,j))/dy;
     elseif XHy_i==inf
        Limiter1y_i=2*(c1(i,j+1)-c1(i,j))/dy;
     elseif XHy_i==-inf
        Limiter1y_i=1*(c1(i,j+1)-c1(i,j))/dy;
     elseif XHy_i<0
        Limiter1y_i=(XHy_i*XHy_i+XHy_i)/(1+XHy_i*XHy_i)*(c1(i,j+1)-c1(i,j))/dy;
     elseif 0<=XHy_i && XHy_i<=1/2
        Limiter1y_i=(XHy_i*XHy_i*XHy_i-2*XHy_i*XHy_i+2*XHy_i)*(c1(i,j+1)-c1(i,j))/dy;
     elseif XHy_i>1/2 && XHy_i<=2
         Limiter1y_i=(3/4*XHy_i+1/4)*(c1(i,j+1)-c1(i,j))/dy;
     elseif XHy_i>2
         Limiter1y_i=(2*XHy_i*XHy_i-2*XHy_i-9/4)/(XHy_i*XHy_i-XHy_i-1)*(c1(i,j+1)-c1(i,j))/dy;
     end


     if isnan(XHx_i_min_1)
        Limiter1x_i_min_1=2*(c1(i,j)-c1(i-1,j))/dx;%%%sigama(i-1)
     elseif XHx_i_min_1==inf
        Limiter1x_i_min_1=2*(c1(i,j)-c1(i-1,j))/dx;
     elseif XHx_i_min_1==-inf
        Limiter1x_i_min_1=1*(c1(i,j)-c1(i-1,j))/dx;
     elseif XHx_i_min_1<0
        Limiter1x_i_min_1=(XHx_i_min_1*XHx_i_min_1+XHx_i_min_1)/(1+XHx_i_min_1*XHx_i_min_1)*(c1(i,j)-c1(i-1,j))/dx;
     elseif 0<=XHx_i_min_1 && XHx_i_min_1<=1/2
        Limiter1x_i_min_1=(XHx_i_min_1*XHx_i_min_1*XHx_i_min_1-2*XHx_i_min_1*XHx_i_min_1+2*XHx_i_min_1)*(c1(i,j)-c1(i-1,j))/dx;
     elseif XHx_i_min_1>1/2 && XHx_i_min_1<=2
        Limiter1x_i_min_1=(3/4*XHx_i_min_1+1/4)*(c1(i,j)-c1(i-1,j))/dx;
     elseif XHx_i_min_1>2
        Limiter1x_i_min_1=(2*XHx_i_min_1*XHx_i_min_1-2*XHx_i_min_1-9/4)/(XHx_i_min_1*XHx_i_min_1-XHx_i_min_1-1)*(c1(i,j)-c1(i-1,j))/dx;%%%sigama(i-1)
     end


     if isnan(XHy_i_min_1)
        Limiter1y_i_min_1=2*(c1(i,j)-c1(i,j-1))/dy;%%%sigama(i-1)
     elseif XHy_i_min_1==inf
        Limiter1y_i_min_1=2*(c1(i,j)-c1(i,j-1))/dy;
     elseif XHy_i_min_1==-inf
        Limiter1y_i_min_1=1*(c1(i,j)-c1(i,j-1))/dy;
     elseif XHy_i_min_1<0
        Limiter1y_i_min_1=(XHy_i_min_1*XHy_i_min_1+XHy_i_min_1)/(1+XHy_i_min_1*XHy_i_min_1)*(c1(i,j)-c1(i,j-1))/dy;
     elseif 0<=XHy_i_min_1 && XHy_i_min_1<=1/2
        Limiter1y_i_min_1=(XHy_i_min_1*XHy_i_min_1*XHy_i_min_1-2*XHy_i_min_1*XHy_i_min_1+2*XHy_i_min_1)*(c1(i,j)-c1(i,j-1))/dy;
     elseif XHy_i_min_1>1/2 && XHy_i_min_1<=2
        Limiter1y_i_min_1=(3/4*XHy_i_min_1+1/4)*(c1(i,j)-c1(i,j-1))/dy;
     elseif XHy_i_min_1>2
        Limiter1y_i_min_1=(2*XHy_i_min_1*XHy_i_min_1-2*XHy_i_min_1-9/4)/(XHy_i_min_1*XHy_i_min_1-XHy_i_min_1-1)*(c1(i,j)-c1(i,j-1))/dy;%%%sigama(i-1)
     end


     if isnan(XHx_i_plu_1)
        Limiter1x_i_plu_1=2*(c1(i+2,j)-c1(i+1,j))/dx;%%%sigama(i+1)
     elseif XHx_i_plu_1==inf
        Limiter1x_i_plu_1=2*(c1(i+2,j)-c1(i+1,j))/dx;
     elseif XHx_i_plu_1==-inf
        Limiter1x_i_plu_1=1*(c1(i+2,j)-c1(i+1,j))/dx;
     elseif XHx_i_plu_1<0
        Limiter1x_i_plu_1=(XHx_i_plu_1*XHx_i_plu_1+XHx_i_plu_1)/(1+XHx_i_plu_1*XHx_i_plu_1)*(c1(i+2,j)-c1(i+1,j))/dx;
     elseif 0<=XHx_i_plu_1 && XHx_i_plu_1<=1/2
        Limiter1x_i_plu_1=(XHx_i_plu_1*XHx_i_plu_1*XHx_i_plu_1-2*XHx_i_plu_1*XHx_i_plu_1+2*XHx_i_plu_1)*(c1(i+2,j)-c1(i+1,j))/dx;
     elseif XHx_i_plu_1>1/2 && XHx_i_plu_1<=2
        Limiter1x_i_plu_1=(3/4*XHx_i_plu_1+1/4)*(c1(i+2,j)-c1(i+1,j))/dx;
     elseif XHx_i_plu_1>2
        Limiter1x_i_plu_1=(2*XHx_i_plu_1*XHx_i_plu_1-2*XHx_i_plu_1-9/4)/(XHx_i_plu_1*XHx_i_plu_1-XHx_i_plu_1-1)*(c1(i+2,j)-c1(i+1,j))/dx;%%%sigama(i+1)
     end


     if isnan(XHy_i_plu_1)
        Limiter1y_i_plu_1=2*(c1(i,j+2)-c1(i,j+1))/dy;%%%sigama(i+1)
     elseif XHy_i_plu_1==inf
        Limiter1y_i_plu_1=2*(c1(i,j+2)-c1(i,j+1))/dy;
     elseif XHy_i_plu_1==-inf
        Limiter1y_i_plu_1=1*(c1(i,j+2)-c1(i,j+1))/dy;
     elseif XHy_i_plu_1<0
        Limiter1y_i_plu_1=(XHy_i_plu_1*XHy_i_plu_1+XHy_i_plu_1)/(1+XHy_i_plu_1*XHy_i_plu_1)*(c1(i,j+2)-c1(i,j+1))/dy;
     elseif 0<=XHy_i_plu_1 && XHy_i_plu_1<=1/2
        Limiter1y_i_plu_1=(XHy_i_plu_1*XHy_i_plu_1*XHy_i_plu_1-2*XHy_i_plu_1*XHy_i_plu_1+2*XHy_i_plu_1)*(c1(i,j+2)-c1(i,j+1))/dy;
     elseif XHy_i_plu_1>1/2 && XHy_i_plu_1<=2
        Limiter1y_i_plu_1=(3/4*XHy_i_plu_1+1/4)*(c1(i,j+2)-c1(i,j+1))/dy;
     elseif XHy_i_plu_1>2
        Limiter1y_i_plu_1=(2*XHy_i_plu_1*XHy_i_plu_1-2*XHy_i_plu_1-9/4)/(XHy_i_plu_1*XHy_i_plu_1-XHy_i_plu_1-1)*(c1(i,j+2)-c1(i,j+1))/dy;%%%sigama(i+1)
     end

    % ���ؽ��
    Limiter1x = [Limiter1x_i; Limiter1x_i_min_1; Limiter1x_i_plu_1];
    Limiter1y = [Limiter1y_i; Limiter1y_i_min_1; Limiter1y_i_plu_1];
end
