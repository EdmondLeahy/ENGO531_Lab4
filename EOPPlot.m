clc
clear
close all
EOP = load('ROPs.txt');
EOP = EOP';
worldPoints = load('Intersection_Pair0_1_FINAL_Obs.txt');
camposes=zeros(6,3,size(EOP,1));
displayscalefactor = 60;
f_l = 18.453986620;
PS_disp = 0.008609300/displayscalefactor;
dx=1728.000000000/2*PS_disp;
dy=2592.000000000/2*PS_disp;

    for tnum=1:size(EOP,1)

        T2=EOP(tnum,1:3);

        R2_ito2=Rotation_e2i(EOP(tnum,4),EOP(tnum,5),EOP(tnum,6))';

        X1=R2_ito2*[0;0;0/displayscalefactor]+T2';

        X2=R2_ito2*[0;0;-f_l/displayscalefactor]+T2';

        X3=R2_ito2*[-dx;dy;-f_l/displayscalefactor]+T2';

        X4=R2_ito2*[dx;dy;-f_l/displayscalefactor]+T2';

        X5=R2_ito2*[dx;-dy;-f_l/displayscalefactor]+T2';

        X6=R2_ito2*[-dx;-dy;-f_l/displayscalefactor]+T2';

        camposes(:,:,tnum)=[X1';X2';X3';X4';X5';X6'];

    end
    figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

    % Create axes
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',18,...
        'FontName','times');
    view(axes1,[0 0]);
    hold(axes1,'all');
   % scatter3(worldPoints(:,1),worldPoints(:,2),worldPoints(:,3),'b*');
    for tnum=1:size(EOP,1)
        %plot3(EOP(tnum,1),EOP(tnum,2),EOP(tnum,3),'gs');
        plot3(camposes(1,1,tnum),camposes(1,2,tnum),camposes(1,3,tnum),'ko','MarkerEdgeColor','k',...
                    'MarkerFaceColor','k',...
                    'MarkerSize',10);
        plot3([camposes(1,1,tnum) camposes(3,1,tnum)],[camposes(1,2,tnum) camposes(3,2,tnum)],[camposes(1,3,tnum) camposes(3,3,tnum)],'k-','LineWidth',1)    
        plot3([camposes(1,1,tnum) camposes(4,1,tnum)],[camposes(1,2,tnum) camposes(4,2,tnum)],[camposes(1,3,tnum) camposes(4,3,tnum)],'k-','LineWidth',1) 
        plot3([camposes(1,1,tnum) camposes(5,1,tnum)],[camposes(1,2,tnum) camposes(5,2,tnum)],[camposes(1,3,tnum) camposes(5,3,tnum)],'k-','LineWidth',1)
        plot3([camposes(1,1,tnum) camposes(6,1,tnum)],[camposes(1,2,tnum) camposes(6,2,tnum)],[camposes(1,3,tnum) camposes(6,3,tnum)],'k-','LineWidth',1)    

        plot3([camposes(3,1,tnum) camposes(4,1,tnum)],[camposes(3,2,tnum) camposes(4,2,tnum)],[camposes(3,3,tnum) camposes(4,3,tnum)],'k-','LineWidth',1)   

        plot3([camposes(4,1,tnum) camposes(5,1,tnum)],[camposes(4,2,tnum) camposes(5,2,tnum)],[camposes(4,3,tnum) camposes(5,3,tnum)],'k-','LineWidth',1)  

        plot3([camposes(5,1,tnum) camposes(6,1,tnum)],[camposes(5,2,tnum) camposes(6,2,tnum)],[camposes(5,3,tnum) camposes(6,3,tnum)],'k-','LineWidth',1)  

        plot3([camposes(6,1,tnum) camposes(3,1,tnum)],[camposes(6,2,tnum) camposes(3,2,tnum)],[camposes(6,3,tnum) camposes(3,3,tnum)],'k-','LineWidth',1)

    end
    %waitfor(figure1);
    
   