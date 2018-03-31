function []= Pelican()
clc
clear all;
close all;

I1=.1213;  I2 = .0116; m1=6.5225; r1=0.0983; m2=2.0458; r2=.0229; l1=.26; l2=.26;

a = I1+I2+m1*r1^2+ m2*(l1^2+ r2^2);
b = m2*l1*r2;
d = I2+ m2*r2^2;

x0= [0.33,0.1,0.0,0.0]; % Initial Condition - Format:[theta1,theta2,dtheta1,dtheta2]

tf=2;

%% Solve the closed-loop system nonlinear differential equation (PlanarArmODE) via ode45
%%ode45 solves the differential equation and returns X with respect to T.
global torque
torque=[];
global sum;
global last_t;
sum = [0;0];
last_t = 0;
[T,X] = ode23(@(t,x)planarArmODE(t,x),[0 tf],x0);

%% Plot Data
figure('Name','Positions under PID SetPoint Control');
plot(T, X(:,1),'r-');
hold on
plot(T, X(:,2),'b--');
hold on

% figure('Name','Input_PD control');
% plot(T, torque(1,1:size(T,1)),'-' );
% hold on
% plot(T, torque(2,1:size(T,1)),'r--');

torque=[];

%% Definging Functions

    function dx = planarArmODE(t,x)
        t
        theta_d=[0;0]; % Desired Set-Point Position
        dtheta_d=[0;0]; % Desired velocity (Derivative of theta_d)
        ddtheta_d=[0;0];
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        
        
        global Mmat Cmat 
        Mmat = [a+2*b*cos(x(2)), d+b*cos(x(2));  d+b*cos(x(2)), d];
        Cmat = [-b*sin(x(2))*x(4), -b*sin(x(2))*(x(3)+x(4)); b*sin(x(2))*x(3),0];
        invM = inv(Mmat);
        invMC = invM*Cmat;
                
        tau = PIDControl(theta_d,dtheta_d,theta,dtheta,t);
        g = 9.81;
        G_matrix = [(m1*r1 + m2*l1)*g * sin(theta(1)) + m2 * r2 * g * sin(theta(1) + theta(2));...
            m2 * r2 * g * sin(theta(1) + theta(2))];
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = -invMC* x(3:4) + invM*tau - invM*G_matrix; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
    end


 function tau = PDControl(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        Kp=100*eye(2);
        Kv=100*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e + Kv*de;
 end
 function tau = PIDControl(theta_d,dtheta_d,theta,dtheta,t)
        Kp=[30,0;...
            0,30];
        Kv=[7,0;...
            0,3];
        Ki=[70,0;...
            0,100];
        e=theta_d-theta; % position error
        dt = (t - last_t)
        if 0==0 
            de = [0;0];
        else
            de = (dtheta_d - dtheta)/dt; % velocity error
            sum = sum + e*dt;
        end
        last_t = t;

        tau = Kp*e + Kv*de + Ki*sum;
    end
    
disp('Finish.');

end
