function jmva_fig3(a, d, Ns)
% Generating the data and initializing the parameters
if nargin<1
    a = 20;
end
if nargin<2
    d = 16;
end
if nargin < 3
    Ns = 10000;
end
if true
    C = rand(d);
    C = C*C'+eye(d);
    Data = eg_rand(C,a,Ns);
    initC = eye(size(C));
    %save sampleD
else 
    load sampleD
end
manopt_methods={'LBFGS','CG','TR'}; %'SD',
nal = 0;
if true
    nal=nal+1;
    if a<d/2
        % Proposed Procedure for non-concave case
        disp('Cest for log-contractive case');
        [Cest,res] = eg_contract_estimate(Data,initC,a);
    else
        % Proposed Procedure for concave case
        disp('Cest for log-concave case');
        [Cest,res] = eg_concave_estimate(Data,initC,a);
    end
    if nal == 1
        begt=res.tocs(1);
    end
    xData{nal} = res.tocs-res.tocs(1)+begt;
    yData{nal} = res.fvals;
    options.legend{nal} = 'fixed-point';
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
        [Cest,res] = manopt_estimate(Data,initC,a,manopt_methods{k});
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

if 1   
    if a<d/2
        nal = nal + 1;
        % Kent-Tyler Procedure for log-contractive case
        disp('Cest for log-contractive case');
        [C,res] = eg_contract_tyler_estimate(Data,initC,a);
        if nal == 1
            begt=res.tocs(1);
        end
        xData{nal} = res.tocs-res.tocs(1)+begt;
        yData{nal} = res.fvals;
        options.legend{nal} = 'Kent-Tyler';
    end
end

%save Datas xData yData alpha beta d b initC Data

options.colors = [0 0 .5
    0.5 0 0
	0 .5 0
    0.5 0.5 0
    0 0.5 0.5
    0.5 0 0.5];
options.xlabel = 'log Running time (seconds) ';
options.ylabel = 'log \Phi(S)-\Phi(S_{min}) ';
options.labelLines = 0;
options.logScale = 3;
options.markerSize = 22;
options.lineWidth = 4;
%options.legendLoc = 'SouthWest';
%options.labelLines = 0;
options.gcaFontSize= 20;
options.titleFontSize = 22;
options.labelFontSize = 22;
options.legendFontSize = 18;
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
for k=1:nal
    p{k}=[pre num2str(round(log10(d(k))*100)/100)];
    lenp(k)=length(p{k});
end
maxlen=max(lenp);
for k=1:nal
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
for k=1:nal
    p{k}=[pre num2str(ceil(log10(dy(k))*100)/100)];
    lenp(k)=length(p{k});
end
maxlen=max(lenp);
for k=1:nal
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
ylabel('Log [l(\Sigma_{max})-l(\Sigma)]','Fontsize',18)
set(gca,'Fontsize',18,'LineWidth',1.8)
set(gcf,'Color',[1 1 1])
if false
    namesave='concave';
    export_fig([namesave '.pdf'],gcf);
    saveas(gcf,namesave, 'png');
    saveas(gcf,namesave, 'fig');
end