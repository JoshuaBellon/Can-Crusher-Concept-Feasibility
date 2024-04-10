clear;
close all;
clc;


l1 = 0.75; %(d) crank arm
l2 = 2; %(a) connecting rod
l4 = 0; %(c) offset
slidMax = sqrt((l1+l2)^2 - l4^2); %max sliding displacement
slidMin = sqrt((l2-l4)^2) - l1; %min sliding displacement
theta1 = linspace(0,pi,181); %span of thetas
%theta1 = linspace(0,2*pi,361);
l3 = slidMax - (sqrt((l1+l2)^2 - l4^2) - ( l1.*cos(theta1) + sqrt((l2^2) - ((l1.*sin(theta1) + l4).^2)) )); %distance bet/ slider and crank pivot
x = sqrt(l3.^2 + l4.^2);
theta3 = asin((l1*sin(theta1))/l2); % angle between connecting rod and ground
theta2 = pi - theta1 - theta3; %angle between connnecting rod and crank arm
%l3 = l2*cos(theta3) + l1*cos(theta1); %(b) %slider displacement assuming no offset


h = 62; %The degree where you engage with the clamp
;

F1 = ones(1,length(theta1));
%This loop should give the graph for the most optimal slider engagement location
for n = 1:length(theta1)
  if (l3(n) - slidMin) - 1 < 0.005 %
    fprintf('Angle Index: %3.2f \n', n);
    i = (length(theta1) - n);
    engage = linspace(0,1,i);
    engage = 7000.*engage;
    for v = 1:i
      F1(n+v) = F1(n+v)*engage(v);
    end
    break
  end
end

%Comment out the loop below to test additional engagement points, you may have to alter the linspace order
%in the loop if you want to flip the direction of the crank slider
%{
for n = h:length(theta1)
  if 1 - (l3(h) - l3(n)) < 0.005 % You have to go within 0.00005 error and change the difference to 0.8922 to get the same result as above
    %fprintf('Please: %3.2f \n', l3(h) - l3(n));
    i = (n - h);
    engage = linspace(0,1,i);
    engage = 7000.*engage;
    for v = 1:i
      F1(h+v) = F1(h+v)*engage(v);
    end
    break
  end
end
%}

%Static Torque Calculations
F2 = F1./cos(theta3);
r = l3.*tan(theta3);
Torque = F2.*r;
worstTorq = (7000./cos(theta3)).*r; %absolute worst loading scenario case scenario


%F3 = abs(F1.*cos(theta2 - (pi/2)));

%Print Max Force and Worst Torque
j = find(abs(pi/2 - theta2) < 0.01);
Max_force = max(Torque);
Worst_Torque = max(worstTorq);
fprintf('Max Force: %3.2f \n', Max_force);
fprintf('Worst Force: %3.2f \n', Worst_Torque);


theta1d = theta1*180/pi; %change radians to degrees

%Print the Plots
figure(1)
plot(theta1d, Torque, 'b*', theta1d, worstTorq, 'k-')
title('Torque needed vs crank angle')

figure(2)
plot(x, Torque, 'm*', x, worstTorq, 'k-')
title('Torque needed vs slide distance')
