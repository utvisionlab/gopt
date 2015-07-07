function jmva_fig2a
qq = [2 3 5 20]; % Dimension
colors = {'b', 'r', 'g', 'k'};
for k2 = 1:2
    if k2 == 1
        eidx = 16;
        bidx = 0.25;
        xticks = [0.25 0.5 1 2 4 8 16];
        figure(11);
    else
        bidx = 9.99999e-6;
        eidx = 0.1;
        xticks = [10^-5 10^-4 10^-3 10^-2 10^-1];
        figure(12);
    end
    for k = 1:4
        q = qq(k);
        a = (bidx:0.03:eidx)*q/2;
        d = q/2;
        y = gammaln(d)+d.*log(2)-gammaln(a)-a*log(q)+ ...
            a.*log(a)+(a-d).*(psi(a)+log(q./a))-a+d;
        if k==1
            if k2==1
                semilogx(a/q*2,y,colors{k},'LineWidth',1.4),hold
            else
                loglog(a/q*2,y,colors{k},'LineWidth',1.4),hold
            end
        else
            plot(a/q*2,y,colors{k},'LineWidth',1.4)
        end
    end
    if k2==1
        set(gca, 'Xlim', [bidx,eidx], 'Xtick', xticks  );
    else
        set(gca,'Xlim',[bidx,eidx],'Xtick', xticks, ...
            'Ytick', [1 10 10^2 10^3 10^4 10^5] );
    end
    set(gca, 'Fontsize', 18, 'LineWidth', 1.8)
    set(gcf, 'Color', [1 1 1])
    ylabel('KL-Divergence', 'Fontsize', 18);
    xlabel('2a/q', 'Fontsize', 18)
    h = legend('Dimension 2', 'Dimension 3', 'Dimension 5', 'Dimension 20');
    set(h,'fontsize',16)
end