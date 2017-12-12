clc;
clear;

% Initialization
size = 400; % Number of grids in x-direction
size2 = 400; % Number of grids in y-direction
counter = 0;
T = zeros(size,size2); % Temperature matrix
p = zeros(size,size2); % Phase field matrix

%Parameter values
del_t = 0.0002;  % Time discretization
del_x = 0.03; % Space discretization
del_y = 0.03; % Space discretization
t = 0;  % Total time taken for simulation
K = 1.6;  % Latent heat
eps = 0.01;
tou = 0.0003;
T_eq = 1;
gamma = 10;
alpha = 0.9;
a = 0.01;

% Initialization of boundaries
for k = 1:10
    p(k,:) = 1; % Makes a portion of the domain on one side solid
end

T(:) = 0;

T_new = zeros(size,size2);
p_new = p;
time = 10;
time_steps = time/del_t;

counter = 0;
l = 0;

%%% Solving for each grid iteratively using explicit scheme %%%
while (l < time_steps)
    for i = 2:size-1
        for j = 2:size2-1
            m = alpha/pi*atan(gamma*(T_eq-T(i,j)));
            if (p(i,j)<=0.5 && p(i,j)>=0)
                p_new(i,j) = p(i,j) + eps^2*del_t/(tou*del_x^2)*(p(i+1,j) + p(i-1,j) + p(i,j+1) +p(i,j-1) - 4*p(i,j)) + del_t/tou*p(i,j)*(1-p(i,j))*(p(i,j)+m-0.5) + a*p(i,j)*(1-p(i,j))*(rand(1)-0.5);
            else
                p_new(i,j) = p(i,j) + eps^2*del_t/(tou*del_x^2)*(p(i+1,j) + p(i-1,j) + p(i,j+1) +p(i,j-1) - 4*p(i,j)) + del_t/tou*p(i,j)*(1-p(i,j))*(p(i,j)+m-0.5);
            end                
            T_new(i,j) = T(i,j) + del_t/del_x^2*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1) -4*T(i,j)) + K*(p_new(i,j)-p(i,j));
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
    T_new(size,:) = T_new(size-1,:) - 100;
    T = T_new;
    p = p_new;
    t = t +del_t;
    l = l + 1;
    counter = counter+1;
    % Plotting at required time intervals
    if (mod(counter,100)==0)
        contourf(p);
        pause(1);
        str1 = ['iso',num2str(counter)];
        print(str1,'-dpng');
    end
    if (mod(counter,1000)==0)
        disp(t);
    end

end

% contourf(p);