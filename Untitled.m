clc
clear


img_ind = load('Inliers_01.txt');

counter_1 = 0;
counter_0 = 0;
for i = 1:size(img_ind)
   
    if img_ind(i) == 1;
        counter_1 = counter_1 + 1;
    end
    if img_ind(i) == 0;
        counter_0 = counter_0 + 1;
    end
end

%% 
close all
clear
int1_xyz = load('Int1_final.txt');
scatter3(int1_xyz(:,1),int1_xyz(:,2),int1_xyz(:,3))

%%

clc
clear
close all

allties = load('AllTies_sparsesift.txt');
allties(:,4) = -allties(:,4);
dlmwrite('AllTies_sparsesift_right.txt',allties,'delimiter','\t','precision',7)
