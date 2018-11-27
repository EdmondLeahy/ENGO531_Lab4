
obs1 = load('obs_int_01.txt');
obs2 = load('obs_int_02.txt');
obs3 = load('obs_int_12.txt');

figure()
hold on
scatter3(obs1(:,1), obs1(:,2), obs1(:,3));
scatter3(obs2(:,1), obs2(:,2), obs2(:,3));
scatter3(obs3(:,1), obs3(:,2), obs3(:,3));

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