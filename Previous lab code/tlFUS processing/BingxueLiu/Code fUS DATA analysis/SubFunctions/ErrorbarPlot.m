%%
[l_data, nbar]=size(data);
data_mean=mean(data,1);
data_std=std(data,1,1);
data_sem=sqrt(std(data,1,1)/sqrt(l_data));

data_err=data_std;

xbar=[-floor(nbar/2):1:floor(nbar/2)];
xbar=[1 1.5 2.5 3];
figure
for ibar=1:nbar
    bar(xbar(ibar),data_mean(ibar),0.4);
    colormap(cool)
    hold on
    h1=errorbar(xbar(ibar),data_mean(ibar),data_err(ibar),'LineStyle','none','Linewidth',1);
    hold on
end
box off
axis tight