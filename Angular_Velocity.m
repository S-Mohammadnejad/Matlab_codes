function [v_t11,v_t12,v_t21,v_t22]=Angular_Velocity(r1,v1,r2,v2)
%Angular_velocity: This function decompose the velocity vector
%in order to find the components which results in w1 and w2
%-------------------------------------------------------------
r12=r1-r2;
v_n1=(dot(v1,r12)/dot(r12,r12)).*r12;
v_t1=v1-v_n1;
v_n2=(dot(v2,r12)/dot(r12,r12)).*r12;
v_t2=v2-v_n2;
while 1
    x1=randn(1,3);
    c11=cross(x1,v_n1);
    c12=dot(x1,v_t1);
    if (c11~=0) & (c12~=0)
        break
    end
end
    
n1=cross(x1,v_n1)/norm(cross(x1,v_n1));
v_t11=(dot(v_t1,n1)/dot(n1,n1)).*n1;
v_t12=v_t1-v_t11;
% while 1
%     x2=randn(1,3);
%     c21=cross(x2,v_n2);
%     c22=dot(x2,v_t2);
%     if (c21~=0) & (c22~=0)
%         break
%     end
% end
n2=cross(x1,v_n2)/norm(cross(x1,v_n2));
v_t21=(dot(v_t2,n2)/dot(n2,n2)).*n2;
v_t22=v_t2-v_t21;
