x = [1,2,2];
y = [2,5,4];
P1 = [x;y]';

x = [0,4,2];
y = [0,4,-6];
P2 = [x;y]';

flag = triangle_intersection(P1, P2);

for i = 1:3
    scatter(P1(i,1),P1(i,2),'r', 'd');
    hold on
    scatter(P2(i,1),P2(i,2),'b');
end
hold off;
grid on

