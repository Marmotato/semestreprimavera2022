z = 0:0.1:20;

J = zeros(5,201);
for i = 0:4
    J(i+1,:) = besselj(i,z);
end

plot(z,J)
grid on
legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
title('Bessel Functions of the First Kind for $\nu \in [0, 4]$','interpreter','latex')
xlabel('z','interpreter','latex')
ylabel('$J_\nu(z)$','interpreter','latex')