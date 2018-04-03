clc
clear all;
close all;
syms t b1 b2 c1 c2 w1 w2


syms q1 q2 q1d q2d
question3_function(t, b1, b2, c1, c2, w1, w2,q1, q2, q1d, q2d)

function []= question3_function(t, b1, b2, c1, c2, w1, w2,q1, q2, q1d, q2d)
%clc
%clear all;
close all;

I1=.1213;  I2 = .0116; m1=6.5225; r1=0.0983; m2=2.0458; r2=.0229; l1=.26; l2=.26;

a = I1+I2+m1*r1^2+ m2*(l1^2+ r2^2);
b = m2*l1*r2;
d = I2+ m2*r2^2;

global Mmat_symb Cmat_symb G_matrix_symb

Mmat_symb = [a+2*b*cos(q1), d+b*cos(q1);  d+b*cos(q1), d];
Cmat_symb = [-b*sin(q1)*q2d, -b*sin(q1)*(q1d+q2d); b*sin(q2)*q1d,0];
g = 9.81;
G_matrix_symb = [(m1*r1 + m2*l1)*g * sin(q1) + m2 * r2 * g * sin(q1 + q2);...
    m2 * r2 * g * sin(q1 + q2)];




x0= [0,0,0.0,0.0]; % Initial Condition - Format:[theta1,theta2,dtheta1,dtheta2]

tf=2;

%% Solve the closed-loop system nonlinear differential equation (PlanarArmODE) via ode45
%%ode45 solves the differential equation and returns X with respect to T.
global torque
torque=[];
global sum;
global last_t;



sum = [0;0];
last_t = 0;
global position_equation
global velocity_equation
global acceleration_equation

e_3t =  exp(-2*(t^3));


position_symb = [...
    b1 * (1-e_3t) + c1 * (1-e_3t) * sin(w1*t);...
    b2 * (1-e_3t) + c2 * (1-e_3t) * sin(w2*t);]

velocity_symb = simplify(diff(position_symb,t))
acceleration_symb = simplify(diff(velocity_symb,t))

position_equation = subs(position_symb,...
    [b1,b2,c1,c2,w1,w2],...
    [pi/4,pi/3,pi/9,pi/6,38.7,118.3]);

velocity_equation = subs(velocity_symb,...
    [b1,b2,c1,c2,w1,w2],...
    [pi/4,pi/3,pi/9,pi/6,38.7,118.3]);

acceleration_equation = subs(acceleration_symb,...
    [b1,b2,c1,c2,w1,w2],...
    [pi/4,pi/3,pi/9,pi/6,38.7,118.3]);
global theta_goals velocity_goals acceleration_goals


[T,X] = ode45(@(t,x)planarArmODE(t,x),[0 tf],x0);

shape = size(T);
theta_goals = ones(shape(1),2);
velocity_goals = ones(shape (1),2);
acceleration_goals = zeros(shape(1),2);
size(theta_goals);

for number = 1:size(1)
    [pos,vel,acc] = Reference_trajectoris(T(number));
    theta_goals(number,:) = [pos(1);pos(2)];
    velocity_goals(number,:) = [vel(1);vel(2)];
    acceleration_goals(number,:) = [acc(1);acc(2)];
end
size(theta_goals);
size(X(:,1:2));
size(T);
error = theta_goals - X(:,1:2);
%% Plot Data
figure('Name','Error ')
plot(T, error(:,1),'r-');
hold on

plot(T, error(:,2),'b--');
hold on


%% Definging Functions

    function dx = planarArmODE(t,x)
        [theta_d,dtheta_d,ddtheta_d] = Reference_trajectoris(t);
        
        theta= x(1:2,1);
        dtheta= x(3:4,1); 
        Mmat = subs(Mmat_symb,[q1;q2;q1d;q2d],x(:,1));
        Cmat = subs(Cmat_symb,[q1;q2;q1d;q2d],x(:,1));
        G_matrix = subs(G_matrix_symb,[q1;q2;q1d;q2d],x(:,1));
        
        Mmat_D = subs(Mmat_symb,[q1;q2;q1d;q2d],[theta_d;dtheta_d]);
        Cmat_D = subs(Cmat_symb,[q1;q2;q1d;q2d],[theta_d;dtheta_d]);
        G_matrix_D = subs(G_matrix_symb,[q1;q2;q1d;q2d],[theta_d;dtheta_d]);

        invMC = Mmat\Cmat;
                
        tau = PDControl(theta_d,dtheta_d,theta,dtheta,t)...
            + Mmat_D*ddtheta_d + Cmat_D * dtheta_d + G_matrix_D;

        %G_matrix = [0;0]; %for when you might as well remove G
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = -invMC* x(3:4) + (Mmat\tau + Mmat\G_matrix) - Mmat\G_matrix; 
    end


 function tau = PDControl(theta_d,dtheta_d,theta,dtheta,t)
        Kp=[200,0;...
           0,150];
        Kv=[3,0;...
          0,3];
        t
        e=theta_d-theta; % position error

        de = dtheta_d - dtheta; % velocity error

        tau = Kp*e + Kv*de;
    end

 function [position,velocity,acceleration] = Reference_trajectoris(time)

 position = subs(position_equation,t,time);
 velocity = subs(velocity_equation,t,time);
 acceleration = subs(acceleration_equation,t,time);
    
 end

    
disp('Finish.');

end
