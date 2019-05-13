%Solve ODEs using Runge Kulta method
%Kara Meyer
%5-12-19

%Using the tic toc command I found that the RK4 method takes almost twice
%as long as the RK2 method but the RK4 method is way more accurate.

%Global varables
omega0=.1;
x0=0; %meters
v0=1; %m/s
t0=0; %s
tf=500; %s
u0=[x0 v0];

%Call functions
tic;
 [t1,u1]=rungek2(@f, t0, tf, 1, u0);
 x1=u1(:,1);
 [t2,u2]=rungek2(@f, t0, tf, 5, u0);
 x2=u2(:,1);
 [t3,u3]=rungek2(@f, t0, tf, 10, u0);
 x3=u3(:,1);
 [t4,u4]=rungek2(@f, t0, tf, 50, u0);
 x4=u4(:,1);
toc;
  
tic;
 [t5,u5]=rungek4(@f, t0, tf, 1, u0);
 x5=u5(:,1);
 [t6,u6]=rungek4(@f, t0, tf, 5, u0);
 x6=u6(:,1);
 [t7,u7]=rungek4(@f, t0, tf, 10, u0);
 x7=u7(:,1);
 [t8,u8]=rungek4(@f, t0, tf, 50, u0);
 x8=u8(:,1);
toc;

%exact solution
texact=t0:.1:tf; %Array of times
xexact=(1/omega0)*sin(omega0*texact);

%Plot the results
figure();
subplot(2,2,1);
plot(texact,xexact,'.',t1,x1,'.',t5,x5,'.');
ylabel('Position (m)');
xlabel('Time (s)');
title('Position over Time - Stepsize=1');

subplot(2,2,2);
plot(texact,xexact,'.',t2,x2,'.',t6,x6,'.');
ylabel('Position (m)');
xlabel('Time (s)');
title('Position over Time - Stepsize=5');

subplot(2,2,3);
plot(texact,xexact,'.',t3,x3,'.',t7,x7,'.');
ylabel('Position (m)');
xlabel('Time (s)');
title('Position over Time - Stepsize=10');
ylim([-50 50]);

subplot(2,2,4);
plot(texact,xexact,'.',t4,x4,'.',t8,x8,'.');
ylabel('Position (m)');
xlabel('Time (s)');
title('Position over Time - Stepsize=50');
ylim([-100 100]);

%Define your ODE
function dudt=f(t,u)

    omega0=.1;

    dudt=zeros(2,1);
    
    dudt(1)=u(2); 
    dudt(2)=-omega0^2*u(1);

end

function [t,u]=rungek2(f,t0, tf, h, initvals)

    %Calculate the number of loops we take
    n=(tf-t0)/h;
    
    %Create an array of zeros to put your values into
    u=zeros(n,2);

    %Set some arrays for time, x, and v values
    u(1,1)=initvals(1);
    u(1,2)=initvals(2);
    t(1)=t0;

    %Loop that will fill in the rest of the v, x, and t array
    for i=2:n
        k1=h*f(t(i-1),u(i-1,:));
        k2=h*f(t(i-1)+.5*h,u(i-1,:)+.5*k1');

        t(i)=t(i-1)+h;
        u(i,:)=u(i-1,:)+k2';
    end

end

function [t,u]=rungek4(f,t0, tf, h, initvals)

    %Calculate the number of loops we take
    n=(tf-t0)/h;
    
    %Create an array of zeros to put your values into
    u=zeros(n,2);

    %Set some arrays for time, x, and v values
    u(1,1)=initvals(1);
    u(1,2)=initvals(2);
    t(1)=t0;

    %Loop that will fill in the rest of the v, x, and t array
    for i=2:n
        k1=h*f(t(i-1),u(i-1,:));
        k2=h*f(t(i-1)+.5*h,u(i-1,:)+.5*k1');
        k3=h*f(t(i-1)+.5*h,u(i-1,:)+.5*k2');
        k4=h*f(t(i-1)+h,u(i-1,:)+k3');

        t(i)=t(i-1)+h;
        u(i,:)=u(i-1,:)+(1/6)*(k1'+2*k2'+2*k3'+k4');
    end

end