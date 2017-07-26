#!/usr/bin/octave

close all
clear all

enable_plot2d=false;
enable_plot3d=true;

% For example with [0, 5], distance to theta=pi/4 should be sqrt(5*5/2)

x=3;
y=2;
z=1; %sqrt(2);

% [x==y] is exactly on the diagonal line pi/4
theta_z=pi/4;

% clockwise rotation around the (non-visible) z-axis
R2z= [ cos(theta_z)  sin(theta_z); 
      -sin(theta_z)  cos(theta_z)];

% counter-clockwise around the (non-visible) z-axis
R2zc=[ cos(theta_z) -sin(theta_z); 
       sin(theta_z)  cos(theta_z)];

p=[x; y];
q=R2z*p;
L=abs(q(2));
printf("Distance: %d\n", L)
f=0;

if (enable_plot2d || enable_plot3d) 
	set(gca(), "defaultlinelinewidth", 4.0);
end

if (enable_plot2d)
	f = f + 1;
	h=figure(f);

	hold on

	% plot the original line
	p0=[0; 0];
	q0=[10; 0];
	q1=R2zc*q0;
	Q=[p0 q1];
	plot(Q(1,:),Q(2,:),'b-');
	
	% plot the original point
	plot(p(1),p(2),'bo');
	
	% plot the projected line
	S=[p0 q0];
	plot(S(1,:),S(2,:),'r-');

	% plot the projected point
	s=R2z*p;
	plot(s(1),s(2),'ro');

	% project the point on the projected line
	sd=s;
	sd(2)=0;
	S1=[s sd];
	plot(S1(1,:),S1(2,:),'r:','linewidth',2);

	% project the original point on the original line
	Sd=[R2zc*sd p];
	plot(Sd(1,:),Sd(2,:),'b:','linewidth',2);

	% calculate the distance
	D=sqrt(sum((Sd(:,1)-Sd(:,2)).^2));
	printf("Calculate distance again: %d\n", D)
	xlabel("x")
	ylabel("y")

	axis("equal")
	hold off
end

% We rotate along another axis, in this case let's choose the y-axis, we rotate from z to x.
theta_y=-pi/6;

% clockwise rotation around the y-axis
R2y= [ cos(theta_y)  sin(theta_y); 
      -sin(theta_y)  cos(theta_y)];

% counter-clockwise around the y-axis
R2yc=[ cos(theta_y) -sin(theta_y); 
       sin(theta_y)  cos(theta_y)];

% Note, in the following R3y depends on R2yc. This might be seen as incorrect. 
% The matrix shifts around such that sin(theta) gets a minus sign -sin(theta) just as with the 2D version and we
% use R2yc rather than R2y.
% See matrices at e.g. Wikipedia, https://en.wikipedia.org/wiki/Rotation_matrix.
% We opted to conform to the Euler angles, which means:
%  R3z rotates from x to y coordinates counter clockwise,
%  R3y rotates from z to x coordinates counter clockwise,
%  R3x rotates from y to z coordinates counter clockwise (not used).
%
% In this way, theta_z>0 and theta_y<0 draws the line through the (+,+,+) positive octant I.
% https://en.wikipedia.org/wiki/Octant_(solid_geometry)

% clockwise matrices
R3x=eye(3);
R3x([2,3],[2,3])=R2z;

R3y=eye(3);
R3y([1,3],[1,3])=R2yc;

R3z=eye(3);
R3z([1,2],[1,2])=R2z;

% counter-clockwise matrices
R3xc=eye(3);
R3xc([2,3],[2,3])=R2zc;

R3yc=eye(3);
R3yc([1,3],[1,3])=R2y;

R3zc=eye(3);
R3zc([1,2],[1,2])=R2zc;

% We assume we have first rotated with theta_z, then with theta_y.
R3c=R3yc*R3zc;

% So we will do here in the reverse order.
% Which is represented by a matrix multiplication.
R3=R3z*R3y;

p=[x; y; z];
q=R3*p;
L=sqrt(q(2).^2+q(3).^2);
printf("Distance in 3D: %d\n", L)

if (enable_plot3d)
	f = f + 1;
	figure(f);

	hold on
	
	% plot the original line (with theta_z and theta_y
	p0=[0; 0; 0];
	q0=[10; 0; 0];
	q2=R3c*q0;
	Q=[p0 q2];
	plot3(Q(1,:),Q(2,:),Q(3,:),'-');

	% plot the original data point
	plot3(p(1),p(2),p(3),'bo');

	% plot the projected line
	S=[p0 q0];
	plot3(S(1,:),S(2,:),S(3,:),'r-');
	
	% plot the projected point
	s=R3*p;
	plot3(s(1),s(2),s(3),'ro');
	
	% project the point on the projected line
	sd=s;
	sd(2)=0;
	sd(3)=0;
	S1=[s sd];
	plot3(S1(1,:),S1(2,:),S1(3,:),'r:','linewidth',2);

	% project the original point on the original line
	Sd=[R3c*sd p];
	plot3(Sd(1,:),Sd(2,:),Sd(3,:),'b:','linewidth',2);

	% calculate the distance
	D=sqrt(sum((Sd(:,1)-Sd(:,2)).^2));
	printf("Calculate distance again: %d\n", D)
	
	% We end up with a line rotate by theta_z and theta_y
	xlabel("x")
	ylabel("y")
	zlabel("z")

	axis("equal")

	az = 80;
	el = 50;
	view(az, el);

	hold off

end
