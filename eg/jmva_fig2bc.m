function jmva_fig2bc
mlambda = 50;
mtheta = 50;
a = 0.5;
q = 2;
Z2 = zeros(mlambda+1,mtheta);
Z1 = zeros(mlambda+1,mtheta);
for k1 = 1:mlambda+1
    lambda = (k1-1) /mlambda *0.9 + 0.1;
    if mod(k1,10) == 0
        fprintf('.');
    end
    for k2 = 1:mtheta+1
        theta = (k2-1)*pi/2 / mtheta;
        covar = [lambda 0; 0 1];
        rot = [cos(theta) sin(theta);-sin(theta) cos(theta)];
        d = eig(covar, rot'*covar*rot, 'chol');
        Z2(k1,k2) = eg_kl(a,covar,a,rot'*covar*rot);
        Z1(k1,k2) = 1/q * sum(d) - 1;
        Z2(k1,k2) = Z2(k1,k2) - a * Z1(k1,k2);
        Z2(k1,k2) = Z2(k1,k2) / (q/2 -a);
    end
end

for k=1:2
    dist1 = 1/mlambda *0.9;
    dist2 = 1/mtheta *pi/2;
    [X,Y] = meshgrid(1./(0.1:dist1:1),(0:dist2:pi/2)/pi);
    if k==1
        figure(21),[c,h] = contour(Y,X,Z2', 'LineWidth', 1.2);
        set(h,'levellist',[0.02 0.1 0.25 0.4 0.55 0.7 0.8 0.9 1 1.1],...
            'textlist',[0.02 0.1 0.25 0.4 0.55 0.7 0.8 0.9 1 1.1])
    else
        figure(22),[c,h] = contour(Y,X,Z1', 'LineWidth', 1.2);
        set(h,'levellist',[0.01 0.1 0.3 0.6 1 1.5 2 2.5 3 3.4 3.8 4],...
            'textlist',[0.01 0.1 0.3 0.6 1 1.5 2 2.5 3 3.4 3.8 4])
    end
    clabel(c,h,'fontsize',16,'color','k')
    set(gca, 'Fontsize', 18, 'LineWidth', 1.8)
    set(gcf, 'Color', [1 1 1])
    xlabel('Rotation degree', 'Fontsize', 18);
    ylabel('\lambda_{max} / \lambda_{min}', 'Fontsize', 18);  
end
%figure,imagesc(Z1),colorbar