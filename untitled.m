
my_pal = [0.4660 0.6740 0.1880 ;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.6350 0.0780 0.1840];


figure

x = [10 9.5 10.2 10]
y = [10 10.2 9.8, 10]
z= [10 10.2 9.8, 9.6]

for i = 1:numel(x)

    scatter3(x(i),y(i),z(i), "MarkerFaceColor", my_pal(i,:), "MarkerEdgeColor", 'k')
    hold on
end

xlim([9, 11])
ylim([9, 11])
zlim([9, 11])

xlabel('X')
ylabel('Y')
zlabel('Z')

grid on