


%%

close all

img_1 = imread('IM_01.jpg');


minx = 0;
miny = 0;
maxx =2592;
maxy =1728;

imagesc([minx maxx], [miny maxy], img_1)

hold on
grid on

obs1 = load('Inliers_Pair0_0.txt');
scatter(obs1(:,1), obs1(:,2), '*')
%obs1 = load('Inliers_Pair0_1.txt');
%scatter(obs1(:,1), obs1(:,2) , '*')

title('2D Points from Ransac Inliers')
xlabel('X (m)')
ylabel('Y (m)')


%%

clc
close all
obs1 = load('Intersection_Pair0_1_FINAL_Obs.txt');
obs2 = load('Intersection_Pair0_2_FINAL_Obs.txt');
obs3 = load('Intersection_Pair1_2_FINAL_Obs.txt');
ROPs = load('ROPs.txt');

figure()
hold on
grid on

scatter3(ROPs(1,1),ROPs(2,1),ROPs(3,1), '*')
scatter3(ROPs(1,2),ROPs(2,2),ROPs(3,2), '*')
scatter3(ROPs(1,3),ROPs(2,2),ROPs(3,3), '*')

scatter3(obs1(:,1), obs1(:,2), obs1(:,3), '*');
scatter3(obs2(:,1), obs2(:,2), obs2(:,3), '*');
scatter3(obs3(:,1), obs3(:,2), obs3(:,3), '*');

title('3D Intersected Points from Three images')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

fprintf('Please click in order: Green board bottom topleft, topright, other left, other right');

[x,y] = ginput(6);

scale = ( dist(x(1),y(1))/0.855 + dist(x(2),y(2))/0.604 + dist(x(3),y(3))/0.351)/3;

obs1(:,1) = obs1(:,1) * (1/scale);
obs1(:,2) = obs1(:,2) * (1/scale);
obs1(:,3) = obs1(:,3) * (1/scale);

obs2(:,1) = obs2(:,1) * (1/scale);
obs2(:,2) = obs2(:,2) * (1/scale);
obs2(:,3) = obs2(:,3) * (1/scale);


figure()
hold on
scatter3(obs1(:,1), obs1(:,2), obs1(:,3));
scatter3(obs2(:,1), obs2(:,2), obs2(:,3));
scatter3(obs3(:,1), obs3(:,2), obs3(:,3));

dist = sqrt(x^2 + y^2);