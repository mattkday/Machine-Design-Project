%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kaleb Coleman, Matthew Day, Scott Meyers, Andrew Peterson
% October 7, 2021
% ME4133
% Project Two
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preliminaries
clear; clc; close all;

%Inital guess for Position Analysis
%Theta 3 is equal to Theta 5
R1=3.52; 
t1=(90*(pi/180));
R2=.82; 
R3=4;
t3=(-110*(pi/180));
R4=1;
t4=0;
R5=2;
t5=(-101*(pi/180));
R6=1.85; 
t6=(63*(pi/180));

%Import .txt data file as a single matrix
File_Data = readmatrix('data');

%Define matrix values to variables to sync two theta and calculate % error
%for data validation
theta2_Pre = File_Data(:,1');
theta3 = File_Data(:,2) ;
Provided_R3 =File_Data(:,3);
Provided_R4 =File_Data(:,4);
Provided_R5 =File_Data(:,5);

%Position Analysis -------------------------------------------------------

%Container for position solutions
M = [];

%Iterators
I = 0;
r = 0;

%Theta 2 setup
theta2 = theta2_Pre';

%Theta 2 iteration from 0 to 2pi
for t2 = theta2; 
    r = r+1;

    %Newton Raphson method
    while (.0005 < R1*exp(1i*(pi/2))+R4*exp(1i*0)+R3*exp(1i*t3)-R2*exp(1i*t2)) | (.0005 < R2*exp(1i*t2)-R5*exp(1i*t3)-R6*exp(1i*t6))
        %^^^^ check if current estimate is close enough using VLE's
        I = I+1;

        % A position matrix
        A = [-R3*sin(t3) cos(t3) 1 0;
          R3*cos(t3) sin(t3) 0 0;
            R5*sin(t3) 0 0 -cos(t3);
            -R5*cos(t3) 0 0 -sin(t3)];
        %b position matrix
        b= -[0+R4+R3*cos(t3)-R2*cos(t2);
            R1+0+R3*sin(t3)-R2*sin(t2);
            R2*cos(t2)-R5*cos(t3)-R6*cos(t6);
            R2*sin(t2)-R5*sin(t3)-R6*sin(t6)];
        j= A\b; %Jacobian

        %add deltas to last estimate to get better estimate
        t3 = t3 + j(1);
        R3 = R3 + j(2);
        R4 = R4 + j(3);
        R5 = R5 + j(4);
    end
    %Results stored in matrix
    M(r,1:5) = [t2,t3,R3,R4,R5];
end
M;
% First Order Coefficients-----------------------------------------------
b1= -[R2*sin(t2);
    -R2*cos(t2);
    -R2*sin(t2);
    R2*cos(t2)];
for iter = 1:size(M,1)
    A = [-M(iter,3)*sin(M(iter,2)) cos(M(iter,2)) 1 0;
        M(iter,3)*cos(M(iter,2)) sin(M(iter,2)) 0 0;
        M(iter,5)*sin(M(iter,2)) 0 0 -cos(M(iter,2));
        -M(iter,5)*cos(M(iter,2)) 0 0 -sin(M(iter,2))];

    b1= -[R2*sin(M(iter,1));
        -R2*cos(M(iter,1));
        -R2*sin(M(iter,1));
        R2*cos(M(iter,1))];
    
    M(iter,6:9) = A\b1;
end
%-------------------------------------------------------------------------

%Data Validation by checking our calculations against the provided Excel
%Sheet data, in the form of percent error
Delta_t2 = (abs(M(:,1) - theta2_Pre) / theta2_Pre) * 100;
Delta_t3 = (abs(M(:,2) - Provided_R3) / Provided_R3) * 100;
Delta_R3 = (abs(M(:,3) - Provided_R3) / Provided_R3) * 100; 
Delta_R4 = (abs(M(:,4) - Provided_R4) / Provided_R4) * 100;
Delta_R5 = (abs(M(:,5) - Provided_R5) / Provided_R5) * 100; 

 %Position anaylsis graph--------------------------------------------

%Define Y-points of graph (t3,r3,r4,r5)
t3= M(:,2)';
R3=M(:,3)';
R4=M(:,4)';
R5=M(:,5)';
%Local minimums and maximums for theta3, R3, R4, and R5 which are also the
%link limits for the unknowns

t3_min = islocalmin(t3);
t3_max = islocalmax(t3);

R3_min = islocalmin(R3);
R3_max = islocalmax(R3);

R4_min = islocalmin(R4);
R4_max = islocalmax(R4);

R5_min = islocalmin(R5);
R5_max = islocalmax(R5);

%Allows multiple data lines on one figure
hold on  

%Plot the Position Analysis data
plot(theta2, t3, 'k', theta2, R3, 'b', theta2, R4, 'm', theta2, R5, 'r')

%Plot the local mins and maxes for data validation (used to compare roots on the First Order Graph)
plot(theta2(t3_min), t3(t3_min), 'k*', theta2(t3_max), t3(t3_max), 'k*', theta2(R3_min), R3(R3_min), 'b*', theta2(R3_max), R3(R3_max), 'b*', theta2(R4_min), R4(R4_min), 'm*', theta2(R4_max), R4(R4_max), 'm*',  theta2(R5_min), R5(R5_min), 'r*', theta2(R5_max), R5(R5_max), 'r*')

%Add legend to make the different lines distinguishable 
legend('theta3 (rad)', 'R3 (in)','R4 (in)','R5 (in)','Local min or max (all lines)')

%Define what the graph is
title('Position Analysis')

%Label the x-axis
xlabel('\theta_2 (rad)')

%Label the y-axis
ylabel('Outputs')

%Add grids to make data easier to read
grid on 
grid minor

%Used to convert the numerical radian to the simplified version original
%Sets the increment of values
set(gca,'XTick',0:pi/4:2*pi) ;

%Defines what the major x-axis grid line should be called
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
%---------------------------------------------------------------------------
%Displaying link limits 
fprintf('Theta 3 Lower Limit %s\n', (min(M(:,2))))
fprintf('Theta 3 Upper Limit %s\n', (max(M(:,2))))

fprintf('Length 3 Lower Limit %s\n',(min(M(:,3))))
fprintf('Length 3 Upper Limit %s\n',(max(M(:,3))))

fprintf('Length 4 Lower Limit %s\n',(min(M(:,4))))
fprintf('Length 4 Upper Limit %s\n',(max(M(:,4))))

fprintf('Length 5 Lower Limit %s\n',(min(M(:,5))))
fprintf('Length 5 Upper Limit %s\n',(max(M(:,5))))




