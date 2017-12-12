clc;
clear;

%Initialization
size = 200;
size2 = 200;
T = zeros(size,size2); % Temperature profile
p = zeros(size,size2); % Phase field
theta = zeros(size,size2);
grad_px = zeros(size,size2);
grad_py = zeros(size,size2);
eps = zeros(size,size2);
eps_der = zeros(size,size2);
grad_eps2x = zeros(size,size2);
grad_eps2y = zeros(size,size2);

%Parameter values
del_t = 0.0002;
del_x = 0.03;
del_y = 0.03;
K = 1.2; % Dimensionless Latent Heat
tou = 0.0003;
T_eq = 1;
gamma = 10;
alpha = 0.9;
a = 0.01;
J = 4; % Mode of anisotropy
delta = 0.05;
eps_bar = 0.01;


T_new = zeros(size,size2);
p_new = p;
time = 10; % Total time for simulation
time_steps = time/del_t;

% Initializing nucleus
for i = 1:size
    for j = 1:size2
        if ((i-size/2)^2 + (j-size2/2)^2 < 20)
            p(i,j) = 1;
        end
    end
end


% for k = 1:3
%     p(k,:) = 1;
% end

t = 0;
l = 0;
test = zeros(size,size2);
contourf(p,1);

% Simulation starts
while (l<time_steps)
    for i = 2:size-1
        for j = 2:size2-1
            grad_px(i,j) = (p(i+1,j)-p(i-1,j))/del_x;
            grad_py(i,j) = (p(i,j+1)-p(i,j-1))/del_y;
            
            if (grad_px == 0)
                if (grad_py(i,j) > 0)
                    theta(i,j) = 0.5*pi;
                elseif(grad_py(i,j) < 0)
                    theta(i,j) = -0.5*pi;
                end
            end
            if (grad_px(i,j) > 0)
                if (grad_py(i,j)<0)
                    theta(i,j) = 2*pi + atan(grad_py(i,j)/grad_px(i,j));
                elseif grad_py(i,j)>0
                    theta(i,j) = atan(grad_py(i,j)/grad_px(i,j));
                end
            end
            if (grad_px(i,j) < 0)
                theta(i,j) = pi + atan(grad_py(i,j)/grad_px(i,j));
            end
%                 theta(i,j) = atan(grad_py(i,j)/grad_px(i,j));
%                 test(i,j) = atan(grad_py(i,j)/grad_px(i,j));
            
            eps(i,j) = eps_bar*(1 + delta*cos(J*theta(i,j)));
            
            eps_der(i,j) = -eps_bar*J*delta*sin(J*theta(i,j));
            grad_eps2x(i,j) = (eps(i+1,j)^2 - eps(i-1,j)^2)/del_x;
            grad_eps2y(i,j) = (eps(i,j+1)^2 - eps(i,j-1)^2)/del_y;
        end
    end
    
    for i = 2:size-1
        for j = 2:size2-1
            part1 = -(eps(i+1,j)*eps_der(i+1,j)*grad_py(i+1,j) - eps(i-1,j)*eps_der(i-1,j)*grad_py(i-1,j))/del_x;
            part2 = (eps(i,j+1)*eps_der(i,j+1)*grad_px(i,j+1) - eps(i,j-1)*eps_der(i,j-1)*grad_px(i,j-1))/del_y;
            part3 = grad_eps2x(i,j)*grad_px(i,j) + grad_eps2y(i,j)*grad_py(i,j);
            m = alpha/pi*atan(gamma*(T_eq-T(i,j)));
            
            if (p(i,j)<=0.5 && p(i,j)>=0)
                p_new(i,j) = p(i,j) + (part1 + part2 + part3 + eps(i,j)^2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-4*p(i,j))/(del_x*del_y) + p(i,j)*(1-p(i,j))*(p(i,j)-0.5+m))*del_t/tou + a*p(i,j)*(1-p(i,j))*(rand(1)-0.5);
            else
                p_new(i,j) = p(i,j) + (part1 + part2 + part3 + eps(i,j)^2*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-4*p(i,j))/(del_x*del_y) + p(i,j)*(1-p(i,j))*(p(i,j)-0.5+m))*del_t/tou;
            end
            
            T_new(i,j) = T(i,j) + (T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)-4*T(i,j))*del_t/del_x^2 + K*(p_new(i,j)-p(i,j));
            
        end
    end
    % Giving no flux boundary conditions
    p_new(1,:) = p_new(2,:);
    p_new(size,:) = p_new(size-1,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,size2) = p_new(:,size2-1);
    T_new(:,1) = T_new(:,2);
    T_new(:,size2) = T_new(:,size2-1);
    T_new(1,:) = T_new(2,:);
    T_new(size,:) = T_new(size-1,:);
    T = T_new;
    p = p_new;
    t = t + del_t;
    l = l + 1;
    %Visualization
%     contourf(p,1);
%     pause(0);
    if (mod(l,100) == 0)
        figure(2);
        contourf(p_new,1);
        pause(1);
        str1 = ['aniso',num2str(l)];
        print(str1,'-dpng');
        disp(t);
    end
end


% contourf(p);




