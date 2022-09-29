clc
clear
load gene
load data_all

data1=data_all(1:21,:);
data2=data_all(22:42,:);
data3=data_all(43:63,:);
data4=data_all(64:84,:);
data5=data_all(85:105,:);
data_new=data1+data2+data3+data4+data5;
data=data_new/5;
data=data_all;
gene=gene;

%%原始网络
node=size(data,2);
sample_num=size(data,1);
figure
s=[];
t=[];
m=1;
DAG=zeros(node);
GGG=zeros(size(data,2));
edge_num=0;
for i=1:length(gene)
    if gene(i,3)~=0
        s(m)=gene(i,1);
        t(m)=gene(i,2);
        m=m+1;
        GGG(gene(i,1),gene(i,2))=1;
        edge_num=edge_num+1;
    end
end
GP=digraph(s,t);
plot(GP)


 [data, mu2, sigma3]= standardizeCols(data);


G=zeros(node,node);
dagvalue=zeros(node,node);

data1=data';
for i=1:node
    for j=i+1:node
        kk=data1(i,:);
        hh=data1(j,:);
        pcc=corr(kk',hh');
        MI=-1/2*log(1-pcc^2);
        dagvalue(i,j)=MI;
        MI_ALL(i,j)=MI;
        if MI>=0.05
            G(i,j)=1;
            %G(j,i)=1;
        end
    end
end
%G=ones(16,16,9);
%求方向
result_DAG_PCB=zeros(node,node);



data2=data;
DD=G(:,:);
nSamples = sample_num;    % the number of samples
nEvals = 20000;       % the maximum number of family evaluations
discrete = 0;         % set to 0 for continuous data
interv = 0;           % set to 0 for observational data
rand('state',0);        % generate data randomly
randn('state',0);
clamped=zeros(sample_num,10); %generate data
penalty = log(nSamples)/2;  % weight of free parameter term

penalty = 1;  
[nSamples,numOfVar]=size(data2); % get the size of the data
DAG_PCB=zeros(numOfVar); % used to record the final network
result_DAG_PCB(:,:)= DAGSearch_test(data2,nEvals,0,penalty,discrete,clamped,DD);



figure
s1=[];
t1=[];
N=node;
m=1;
DAG2=result_DAG_PCB(:,:);
DAG3=DAG2;
for i=1:N
    for j=1:N
        if DAG3(i,j)~=0
            s1(m)=i;
            t1(m)=j;
            m=m+1;
        end
    end
end
GP2=digraph(s1,t1);
plot(GP2)
point_num=node;
TP=0;
FP=0;
TN=0;
FN=0;
aaa=0;
bbb=0;
for i=1:point_num
    for j=1:point_num
        if GGG(i,j)==1 && DAG3(i,j)==1
            TP=TP+1;
        elseif GGG(i,j)==0 && DAG3(i,j)==1
            FP=FP+1;
        else 
            aaa=aaa+1;

        end
    end
end
    
fanxiang=aaa+bbb
TP
FN=sum(sum(GGG))-TP
FP
TN=point_num*(point_num-1)-sum(sum(GGG))-FP
Precision = TP/(TP+FP)
TPR = TP/(TP+FN)
FPR = FP/(FP+TN)
Specificity= TN/(FP+TN)
%= 1 - FPR 
Accuracy =(TP+TN)/(TP+TN+FP+FN)
error= (FP+FN)/(TP+TN+FP+FN)

P_XY_Z=TP/(TP+FN);
P_Z=(TP+FN)/(TP+TN+FP+FN);
P_XY=(TP+FP)/(TP+TN+FP+FN);
posterior_probability=(P_XY_Z*P_Z)/P_XY



numberTC=sum(DAG3(:)~=GGG(:))    
edge_num

AAAA(1)=TP;
AAAA(2)=FN ;
AAAA(3)=FP ;
AAAA(4)=TN ;
AAAA(5)=Precision ;
AAAA(6)=TPR ;
AAAA(7)=FPR ;
AAAA(8)=Specificity ;
AAAA(9)=Accuracy ;
AAAA(10)=error ;
AAAA(11)=posterior_probability;