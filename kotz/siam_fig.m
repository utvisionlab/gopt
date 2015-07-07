function siam_fig(a, c, d, Ns, Data, initC)
%
%   siam_fig(a, c, d, Ns)
%
%  This function is for generating comparison plots as in siam paper
%     
%  of Kotz-type distribution (KTD) which has the following density:
%
%  f(X)=Gamma(0.5*n)/(pi^(0.5*n))*c/Gamma(a)*1/b^a*det(C)^-0.5*u^(a*c-n/2)*exp(-U^c/b)
%
%  where  U=X' C^{-1} X,  n  is the dimension.
%  The parameters of the distribution are the nxn Covariance matrix  C 
%  and the shape parameters  a, c .  Ns  is the number of samples
%  Since C is true covariance b is equal to 
%      b=(n*gamma(a)/gamma(a+1/c))^(c)
%  The mean value of the generalized gamma is equal to
%     gamma(a+1/c)/gamma(a)*b^(1/c) 
%  
%  parameter c is the same as parameter beta in the paper
%  parameter a is different than alpha, alpha = a * c
%
if nargin < 3
    d = 16;
end
if nargin < 2
    c = 0.5;
end
if nargin < 1
    a = d/c/4;
end
if nargin < 4
    Ns = 10000;
end

if nargin < 5
    C = rand(d);
    C = C*C'+eye(d);
    Data = kotz_rand(a,c,C,Ns);
end

if nargin < 6
    initC = eye(d);
end

b = (d*gamma(a)/gamma(a+1/c))^c;
func = @(x)2*1/b*(x.^c)+2*(d/2-a.*c).*log(x);
gfunc = @(x)2*c/b*(x.^(c-1))+2*(d/2-a.*c)./x;

manopt_methods={'LBFGS','CG','SD','TR'}; %'SD',
nal = 0;
if true
    nal=nal+1;
    % Proposed Procedure for concave case
    disp('Cest for log-contractive case');
    [Cest,res] = ecd_contract_estimate(Data,initC,func,gfunc);   
    if nal == 1
        begt=res.tocs(1);
    end
    xData{nal} = res.tocs-res.tocs(1)+begt;
    yData{nal} = res.fvals;
    options.legend{nal} = 'FP';
end
if false
    Cest
    C
    Data*Data'/size(Data,2)
end
if true
    for k = 1:numel(manopt_methods)
        nal = nal + 1;
        % Estimating procedure for manopt case
        disp('Cest for manopt case');
        [Cest,res] = manopt_contract_estimate(Data,initC,func,gfunc,manopt_methods{k});
        tocs = zeros(1,numel(res));
        fvals = zeros(1,numel(res));
        for k2 = 1:numel(res)
                tocs(k2) = res(k2).time;
                fvals(k2)= res(k2).cost;
        end
        if nal==1
            begt=tocs(1);
        end
        xData{nal}=tocs-tocs(1)+begt;
        yData{nal}=fvals;
        options.legend{nal}=manopt_methods{k};
    end
end

if true
    nal = nal + 1;
    % Fast Procedure for log-contractive case
    disp('Cest for fast log-contractive case');
    [C,res] = kotz_fast_contract_estimate(Data,initC,a,c);
    if nal == 1
        begt=res.tocs(1);
    end
    xData{nal} = res.tocs-res.tocs(1)+begt;
    yData{nal} = res.fvals;
    options.legend{nal} = 'FP2';
end

if true
    nal=nal+1;
    % Proposed Procedure for concave case
    disp('Cest for line-search log-contractive case');
    [Cest,res] = ecd_ls_contract_estimate(Data,initC,func,gfunc);    
    if nal == 1
        begt=res.tocs(1);
    end
    xData{nal} = res.tocs-res.tocs(1)+begt;
    yData{nal} = res.fvals;
    options.legend{nal} = 'FP-S';
end

%save Datas xData yData alpha beta d b initC Data

options.colors = [0 0 .5
    0.5 0 0
	0 .5 0
    0.5 0.5 0
    0 0.5 0.5
    0.5 0 0.5
    0.33 0.33 0.33
    ];
options.xlabel = 'log Running time (seconds) ';
options.ylabel = 'log \Phi(S)-\Phi(S_{min}) ';
options.labelLines = 0;
options.logScale = 3;
options.markerSize = 22;
options.lineWidth = 4;
%options.legendLoc = 'NorthEast';
%options.labelLines = 0;
options.gcaFontSize= 20;
options.titleFontSize = 22;
options.labelFontSize = 22;
options.legendFontSize = 14;
options.axesLineWidth = 3;
options.gcaTickLength = 0.03;
xmin=Inf;
xmax=-Inf;
ymin=Inf;
ymax=-Inf;
for i=1:length(xData)
    m=min(xData{i}(:));
    xmin=min([xmin,m]);
    m=max(xData{i}(:));
    xmax=max([xmax,m]);
    m=min(yData{i}(:));
    ymin=min([ymin,m]);
    m=max(yData{i}(:));
    ymax=max([ymax,m]);
end
for i=1:length(xData)
    yData{i}=yData{i}-ymin;

end
ymax=ymax-ymin;
ymin=0.00001;
prettyPlot(xData,yData,options);
set(gcf,'Color',[1 1 1])
set(gca,'Xlim',[xmin/1.001 xmax*1.001],'Ylim',[ymin/1.001 ymax*1.001]);
%pre='10^';
pre='';
d=logspace(ceil(log10(xmin)*10)/10, floor(log10(xmax)*10)/10,6);
for k=1:6
    p{k}=[pre num2str(round(log10(d(k))*100)/100)];
    lenp(k)=length(p{k});
end
maxlen=max(lenp);
for k=1:6
    disb=[];
    disa=[];
    for m=1:floor((maxlen-lenp(k))/2)
        disb=[disb ' '];
    end
    for m=1:ceil((maxlen-lenp(k))/2)
        disa=[disa ' '];
    end
    xTickLabel(k,:)=[disb p{k} disa];
end
%set(gca,'Xlim',[xmin/d(3)*d(1) xmax*1.001]);
dy=logspace(ceil(log10(ymin)*10)/10, floor(log10(ymax)*10)/10,6);
for k=1:6
    p{k}=[pre num2str(ceil(log10(dy(k))*100)/100)];
    lenp(k)=length(p{k});
end
maxlen=max(lenp);
for k=1:6
    disb=[];
    disa=[];
    for m=1:floor((maxlen-lenp(k))/2)
        disb=[disb ' '];
    end
    for m=1:ceil((maxlen-lenp(k))/2)
        disa=[disa ' '];
    end
    yTickLabel(k,:)=[disb p{k} disa];
end

set(gca,'Xtick',sort(d),'XTickLabel',xTickLabel,'Ytick',sort(dy),'YTickLabel',yTickLabel);

set(gca,'Fontsize',13,'LineWidth',1.8)
xlabel('Log running time (seconds)','Fontsize',18);
ylabel('Log [l_n(\Sigma_{max})-l_n(\Sigma)]','Fontsize',18)
set(gca,'Fontsize',18,'LineWidth',1.8)
set(gcf,'Color',[1 1 1])
if false
    namesave='siam';
    export_fig([namesave '.pdf'],gcf);
    saveas(gcf,namesave, 'png');
    saveas(gcf,namesave, 'fig');
end