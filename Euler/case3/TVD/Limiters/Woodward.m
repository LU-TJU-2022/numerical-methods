function [Limiter1x, Limiter1y] = Woodward(c1, dx, dy, i, j)
    % Woodward ���������㺯��
    % ����:
    %   c1  - ���������ϵ�Ũ�ȳ�����
    %   dx, dy - ������
    %   i, j - ��ǰ����������
    % ���:
    %   Limiter1x - x ����� Woodward ������ֵ
    %   Limiter1y - y ����� Woodward ������ֵ
    
    XHx_i=(c1(i,j)-c1(i-1,j))/(c1(i+1,j)-c1(i,j));%%%seita(i)
    XHy_i=(c1(i,j)-c1(i,j-1))/(c1(i,j+1)-c1(i,j));%%%seita(i)
    %%%%%%%����limiter����i-1
    XHx_i_min_1=(c1(i-1,j)-c1(i-2,j))/(c1(i,j)-c1(i-1,j));%%%seita(i-1)
    XHy_i_min_1=(c1(i,j-1)-c1(i,j-2))/(c1(i,j)-c1(i,j-1));%%%seita(i-1)
    %%%%%%%%%����limiter����i+1
    XHx_i_plu_1=(c1(i+1,j)-c1(i,j))/(c1(i+2,j)-c1(i+1,j));%%%seita(i+1)
    XHy_i_plu_1=(c1(i,j+1)-c1(i,j))/(c1(i,j+2)-c1(i,j+1));%%%seita(i+1)
    
    Limiter1x_i=max(0,min(min(2,2*XHx_i),0.5*(1+XHx_i)))*(c1(i+1,j)-c1(i,j))/dx;
    Limiter1y_i=max(0,min(min(2,2*XHy_i),0.5*(1+XHy_i)))*(c1(i,j+1)-c1(i,j))/dy;
    Limiter1x_i_min_1=max(0,min(min(2,2*XHx_i_min_1),0.5*(1+XHx_i_min_1)))*(c1(i,j)-c1(i-1,j))/dx;%%%sigama(i-1)
    Limiter1y_i_min_1=max(0,min(min(2,2*XHy_i_min_1),0.5*(1+XHy_i_min_1)))*(c1(i,j)-c1(i,j-1))/dy;%%%sigama(i-1)
    Limiter1x_i_plu_1=max(0,min(min(2,2*XHx_i_plu_1),0.5*(1+XHx_i_plu_1)))*(c1(i+2,j)-c1(i+1,j))/dx;%%%sigama(i+1)
    Limiter1y_i_plu_1=max(0,min(min(2,2*XHy_i_plu_1),0.5*(1+XHy_i_plu_1)))*(c1(i,j+2)-c1(i,j+1))/dy;%%%sigama(i+1)

    % ���ؽ��
    Limiter1x = [Limiter1x_i; Limiter1x_i_min_1; Limiter1x_i_plu_1];
    Limiter1y = [Limiter1y_i; Limiter1y_i_min_1; Limiter1y_i_plu_1];
end
