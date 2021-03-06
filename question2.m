clc
clear all;
close all;
syms t b1 b2 c1 c2 w1 w2
question2_function(t, b1, b2, c1, c2, w1, w2)

function []= question2_function(t, b1, b2, c1, c2, w1, w2)
%clc
%clear all;
close all;

I1=.1213;  I2 = .0116; m1=6.5225; r1=0.0983; m2=2.0458; r2=.0229; l1=.26; l2=.26;

a = I1+I2+m1*r1^2+ m2*(l1^2+ r2^2);
b = m2*l1*r2;
d = I2+ m2*r2^2;






x0= [0,0,0.0,0.0]; % Initial Condition - Format:[theta1,theta2,dtheta1,dtheta2]

tf=10;

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

e_3t =  exp(-2*(t^3))


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
        Mmat = [a+2*b*cos(x(2)), d+b*cos(x(2));  d+b*cos(x(2)), d];
        Cmat = [-b*sin(x(2))*x(4), -b*sin(x(2))*(x(3)+x(4)); b*sin(x(2))*x(3),0];
        invMC = Mmat\Cmat;
                
        control = Computed_Torque(theta_d,dtheta_d,ddtheta_d,theta,dtheta,t);
        g = 9.81;
        G_matrix = [(m1*r1 + m2*l1)*g * sin(x(1)) + m2 * r2 * g * sin(x(1) + x(2));...
            m2 * r2 * g * sin(x(1) + x(2))];
        %G_matrix = [0;0]; %for when you might as well remove G
        
        tau = (control + Cmat * x(3:4) + G_matrix);
        %torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = -invMC* x(3:4) + (tau) - Mmat\G_matrix; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
        %dx(5:6) = theta_d;
    end

function tau = Computed_Torque(theta_d,dtheta_d,ddtheta_d,theta,dtheta,time)
    Kp=[1500,0;...
        0,14000];
    Kv=[77.46,0;...
        0,236.64];

    time
    e=theta_d-theta; % position error
    de = dtheta_d - dtheta; % velocity error

    tau = Kp*e + Kv*de + ddtheta_d;

end

 function [position,velocity,acceleration] = Reference_trajectoris(time)

 position = subs(position_equation,t,time);
 velocity = subs(velocity_equation,t,time);
 acceleration = subs(acceleration_equation,t,time);
    
 end

 function tau = PIDControl(theta_d,dtheta_d,theta,dtheta,t)
        Kp=[30,0;...
           0,30];
        Kv=[7,0;...
          0,3];
        Ki=[70,0;...
            0,100];
        t
        e=theta_d-theta; % position error
        dt = (t - last_t);

        de = dtheta_d - dtheta; % velocity error
        sum = sum + e*dt;

        last_t = t;
        tau = Kp*e + Kv*de + Ki*sum;
    end
    
disp('Finish.');

end
