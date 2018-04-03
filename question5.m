syms q1 q2 q1d q2d
question4_func(q1, q2, q1d, q2d)
function []= question4_func(q1, q2, q1d, q2d)
clc
%clear all;
close all;

I1=1;  I2 = 1; m1=0.5; r1=0.15; m2=0.5; r2=.15; l1=0.3; l2=0.3;

a = I1+I2+m1*r1^2+ m2*(l1^2+ r2^2);
b = m2*l1*r2;
d = I2+ m2*r2^2;

global Mmat_symb Cmat_symb G_matrix_symb


%% Model equations
Mmat_symb = [a+2*b*cos(q1), d+b*cos(q1);  d+b*cos(q1), d];
Cmat_symb = [-b*sin(q1)*q2d, -b*sin(q1)*(q1d+q2d); b*sin(q2)*q1d,0];
g = 9.81;
G_matrix_symb = [(m1*r1 + m2*l1)*g * sin(q1) + m2 * r2 * g * sin(q1 + q2);...
    m2 * r2 * g * sin(q1 + q2)];




%xd = input('X coordinate ');
%yd = input('Y coordinate ');
%% position calculations
posi = [0.300,0.450]
posf = [-0.300,0.450]


[qi1,qi2] = twoDOFIK(posi(1),posi(2),1);
[qf1,qf2] = twoDOFIK(posf(1),posf(2),-1);



x0= [qi1,qi2,0,0]; % Initial Condition - Format:[theta1,theta2,dtheta1,dtheta2]


tf=10;

%% Solve the closed-loop system nonlinear differential equation (PlanarArmODE) via ode45
%%ode45 solves the differential equation and returns X with respect to T.
global torque
torque=[];

[T,X] = ode45(@(t,x)planarArmODE(t,x),[0 tf],x0);

temp = ones(length(X(:,1)));

[x,y] = twoDOFFK([X(:,1),X(:,2)],l1,l2);
poses = [x,y];
linkx = l1*cos(X(:,1));
linky = l1*sin(X(:,1)) ;
links = [linkx,linky];


error_joints = [temp * qf1 - X(:,1),temp * qf1 - X(:,2)];
error_task = sqrt((posf(1)-x).^2 + (posf(2)-y).^2) ;

task_velocity = diff(poses);

%% Plot Data
figure('Name','Model');
plot([0,links(end,1),poses(end,1)],[0,links(end,2),poses(end,2)],'r-');
hold on
plot([0,links(1,1),poses(1,1)],[0,links(1,2),poses(1,2)],'b-');
hold on

figure('Name','Task Error');
plot(T, error_task,'r-');
hold on

figure('Name','Positions red is x');
plot(T, x,'r-');
hold on
plot(T, y,'b-');
hold on

figure('Name','Theta_1 Error');
plot(T, error_joints(:,1),'r-');
hold on

figure('Name','Theta_2 Error');
plot(T, error_joints(:,2),'r-');
hold on

figure('Name','Theta_1 Position');
plot(T, X(:,1),'r-');
hold on

figure('Name','Theta_1 Velocity ');
plot(T, [0;diff(X(:,1))],'r-');
hold on
%}
figure('Name','Theta_2 under PD SetPoint Control');
plot(T, X(:,2),'r--');
hold on
figure('Name','Theta_2 Velocity ');
plot(T, [0;diff(X(:,2))],'r-');
hold on

figure('Name','Input_PD control');
plot(T, torque(1,1:size(T,1)),'-' );
hold on
plot(T, torque(2,1:size(T,1)),'r--');

torque=[];

%% Definging Functions
    function [q1,q2] = twoDOFIK(x,y,direction)
        q2 = direction * acos((x^2 + y^2 -l1^2 - l2^2)/(2*l1*l2));
        B = atan2((l2*sin(q2)),(l1 + l2*cos(q2)));
        Y = atan2(y,x);
        q1 = Y - direction * B;       
    end
    function [x,y] = twoDOFFK(data,l1,l2)
        q1 = data(:,1);
        q2 = data(:,2);
        
        x = l1 * cos(q1) + l2*cos(q1 + q2);
        y = l1 * sin(q1) + l2*sin(q1 + q2);
    end
    function dx = planarArmODE(t,x)
        theta_d=[qf1;-qf2]; % Desired Set-Point Position
        dtheta_d=[0;0]; % Desired velocity (Derivative of theta_d)
        ddtheta_d=[0;0];
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        
        t
        Mmat = subs(Mmat_symb,[q1;q2;q1d;q2d],x(:,1));
        Cmat = subs(Cmat_symb,[q1;q2;q1d;q2d],x(:,1));
        G_matrix = subs(G_matrix_symb,[q1;q2;q1d;q2d],x(:,1));
        
        Mmat_D = subs(Mmat_symb,[q1;q2;q1d;q2d],[theta_d;dtheta_d]);
        Cmat_D = subs(Cmat_symb,[q1;q2;q1d;q2d],[theta_d;dtheta_d]);
        G_matrix_D = subs(G_matrix_symb,[q1;q2;q1d;q2d],[theta_d;dtheta_d]);
        tau = Computed_Torque(theta_d,dtheta_d,ddtheta_d,theta,dtheta,t);
%          tau = PDControl(theta_d,dtheta_d,theta,dtheta,t)...
%              + Mmat_D*ddtheta_d + Cmat_D * dtheta_d + G_matrix_D;
        
            
        
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        %dx(3:4) =   Mmat\(tau-Cmat*x(3:4) - G_matrix);
        dx(3:4) = tau
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
 function tau = PDControl(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        Kp=10*eye(2); 
        Kv=10*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e + Kv*de;
    end
    
disp('Finish.');
end
