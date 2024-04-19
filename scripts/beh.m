%% Beh

low_pass = 0.5; % preprocessed data with two sets of parameters 0.1-30Hz and 0.5-30Hz
savepath = '';
savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save ERPs
%% Behavioral exploration
% 2021-11-30
subj_set = {'ST_04','ST_05','ST_06','ST_09','ST_10','ST_11','ST_12','ST_13','ST_14','ST_15','ST_16','ST_18',...
            'ST_19','ST_20','ST_21','ST_22','ST_23','ST_24','ST_25','ST_26','ST_27','ST_28','ST_29','ST_30','ST_31',...
            'ST_32','ST_33','ST_34','ST_35','ST_36'};
% Difficulty rating
for i=1:length(subj_set)
    filename = [savepath,strcat(subj_set(i),'\beh_',subj_set(i),'_main_coded.xlsx')];
    [~,~,raw] = xlsread(filename{:});
    
    %% check the order make sure no error in coded excel
    cof_f = str2double({raw{2:73,6}}); % confi pic folder
    ret_f = [raw{2:73,11}]; % encoding pic folder
    ret_ind = [raw{2:73,8}]; % ret index

    
    %% Get the beh data
    % enc trial index
    enc_trial_index(i,:) = [raw{2:73,2}];
    
    diff_rating(i,:) = cof_f(enc_trial_index(i,:)); % re-order corresponding to ret order
end

savepath = '';
savename = strcat(savepath,'Beh&index_diffrate_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
%save (savename,'diff_rating','-v7.3');%save ERPs
load(savename)

%% Localizer
savename = strcat(savepath,'Beh&index_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);

%% Number of recall
% number of item
M = mean(num_recall,2,'omitnan');

% n_recall
n_all = zeros(size(num_recall,1),3);
n_all(:,1)=sum(num_recall==0,2)./sum(~isnan(num_recall),2)*100;
n_all(:,2)=sum(num_recall==1,2)./sum(~isnan(num_recall),2)*100;
n_all(:,3)=sum(num_recall==2,2)./sum(~isnan(num_recall),2)*100;
%n_all(13,:)=[];

mean_per = mean(n_all/100,1);
SD_all = std(n_all/100,0,1);

figure(1);clf;
for fig=1
err = std(n_all);
b=bar(mean(n_all,1),0.8);hold on
er = errorbar(mean(n_all,1),err); 
scatter(ones(size(n_all(:,1))).*(1+(rand(size(n_all(:,1)))-0.5)/5),n_all(:,1),'k','filled');
scatter(ones(size(n_all(:,2))).*(2+(rand(size(n_all(:,2)))-0.5)/5),n_all(:,2),'k','filled');
scatter(ones(size(n_all(:,3))).*(3+(rand(size(n_all(:,3)))-0.5)/5),n_all(:,3),'k','filled');hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;
b.FaceColor = 'flat';
b.FaceAlpha = 0.6;
    b.CData(1,:) = [0,0,1];
    b.CData(2,:) = [1,0,0];
xticks([1 2 3])
xticklabels({'No rem','rem 1','rem all'})
ylim([0,80])
ylabel('Percentage of trials')
title('Image Recall')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end 


% n of items recalled Separate by block number
item_blocknumber = zeros(size(num_recall,1),6);
item_blocknumber(:,1) = sum(num_recall(:,1:12),2,'omitnan')./sum(~isnan(num_recall(:,1:12)),2);
item_blocknumber(:,2) = sum(num_recall(:,13:24),2,'omitnan')./sum(~isnan(num_recall(:,13:24)),2);
item_blocknumber(:,3) = sum(num_recall(:,25:36),2,'omitnan')./sum(~isnan(num_recall(:,25:36)),2);
item_blocknumber(:,4) = sum(num_recall(:,36:48),2,'omitnan')./sum(~isnan(num_recall(:,36:48)),2);
item_blocknumber(:,5) = sum(num_recall(:,49:60),2,'omitnan')./sum(~isnan(num_recall(:,49:60)),2);
item_blocknumber(:,6) = sum(num_recall(:,61:72),2,'omitnan')./sum(~isnan(num_recall(:,61:72)),2);

figure(2);clf;
for fig=1
    err = std(item_blocknumber);
    b=bar(mean(item_blocknumber,1),0.8);hold on
    er = errorbar(mean(item_blocknumber,1),err);
    scatter(ones(size(item_blocknumber(:,1))).*(1+(rand(size(item_blocknumber(:,1)))-0.5)/5),item_blocknumber(:,1),'k','filled');
    scatter(ones(size(item_blocknumber(:,2))).*(2+(rand(size(item_blocknumber(:,2)))-0.5)/5),item_blocknumber(:,2),'k','filled');
    scatter(ones(size(item_blocknumber(:,3))).*(3+(rand(size(item_blocknumber(:,3)))-0.5)/5),item_blocknumber(:,3),'k','filled');
    scatter(ones(size(item_blocknumber(:,4))).*(4+(rand(size(item_blocknumber(:,4)))-0.5)/5),item_blocknumber(:,4),'k','filled');
    scatter(ones(size(item_blocknumber(:,5))).*(5+(rand(size(item_blocknumber(:,5)))-0.5)/5),item_blocknumber(:,5),'k','filled');
    scatter(ones(size(item_blocknumber(:,6))).*(6+(rand(size(item_blocknumber(:,6)))-0.5)/5),item_blocknumber(:,6),'k','filled');hold off
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 3;
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
    %b.CData([1,4],:) = [0,0,1;0,0,1];
    %b.CData([2,5],:) = [1,0,0;1,0,0];
    %b.CData([3,6],:) = [0,1,0;0,1,0];
    xticks([1 2 3 4 5 6])
    xticklabels({'B1','B2','B3','B4','B5','B6'})
    ylim([0,2])
    ylabel('Numbers of items')
    title('Image Recall')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end


% n of items recalled Separate by block type
item_blocktype = zeros(size(num_recall,1),6);
for s=1:30
    item_blocktype(s,1) = sum(num_recall(s,block_seq(s,:)==1),'omitnan')/sum(~isnan(num_recall(s,block_seq(s,:)==1)));
    item_blocktype(s,2) = sum(num_recall(s,block_seq(s,:)==2),'omitnan')/sum(~isnan(num_recall(s,block_seq(s,:)==2)));
    item_blocktype(s,3) = sum(num_recall(s,block_seq(s,:)==3),'omitnan')/sum(~isnan(num_recall(s,block_seq(s,:)==3)));
    item_blocktype(s,4) = sum(num_recall(s,block_seq(s,:)==4),'omitnan')/sum(~isnan(num_recall(s,block_seq(s,:)==4)));
    item_blocktype(s,5) = sum(num_recall(s,block_seq(s,:)==5),'omitnan')/sum(~isnan(num_recall(s,block_seq(s,:)==5)));
    item_blocktype(s,6) = sum(num_recall(s,block_seq(s,:)==6),'omitnan')/sum(~isnan(num_recall(s,block_seq(s,:)==6)));

end 


figure(3);clf;
for fig=1
    err = std(item_blocktype);
    b=bar(mean(item_blocktype,1),0.8);hold on
    er = errorbar(mean(item_blocktype,1),err);
    scatter(ones(size(item_blocktype(:,1))).*(1+(rand(size(item_blocktype(:,1)))-0.5)/5),item_blocktype(:,1),'k','filled');
    scatter(ones(size(item_blocktype(:,2))).*(2+(rand(size(item_blocktype(:,2)))-0.5)/5),item_blocktype(:,2),'k','filled');
    scatter(ones(size(item_blocktype(:,3))).*(3+(rand(size(item_blocktype(:,3)))-0.5)/5),item_blocktype(:,3),'k','filled');
    scatter(ones(size(item_blocktype(:,4))).*(4+(rand(size(item_blocktype(:,4)))-0.5)/5),item_blocktype(:,4),'k','filled');
    scatter(ones(size(item_blocktype(:,5))).*(5+(rand(size(item_blocktype(:,5)))-0.5)/5),item_blocktype(:,5),'k','filled');
    scatter(ones(size(item_blocktype(:,6))).*(6+(rand(size(item_blocktype(:,6)))-0.5)/5),item_blocktype(:,6),'k','filled');hold off
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 3;
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
    %b.CData([1,4],:) = [0,0,1;0,0,1];
    %b.CData([2,5],:) = [1,0,0;1,0,0];
    %b.CData([3,6],:) = [0,1,0;0,1,0];
    xticks([1 2 3 4 5 6])
    xticklabels({'F-P-O','P-F-O','O-F-P','F-O-P','P-O-F','O-P-F'})
    ylim([0,2])
    ylabel('Numbers of items')
    title('Image Recall')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end



% Separate by block type
% n_recall
% n of items recalled Separate by block type
rem_typeblock0 = zeros(size(num_recall,1),6);
rem_typeblock1 = zeros(size(num_recall,1),6);
rem_typeblock2 = zeros(size(num_recall,1),6);
for s=1:30
    % rem0
    rem_typeblock0(s,1) = sum(num_recall(s,block_seq(s,:)==1)==0)/sum(~isnan(num_recall(s,block_seq(s,:)==1)));
    rem_typeblock0(s,2) = sum(num_recall(s,block_seq(s,:)==2)==0)/sum(~isnan(num_recall(s,block_seq(s,:)==2)));
    rem_typeblock0(s,3) = sum(num_recall(s,block_seq(s,:)==3)==0)/sum(~isnan(num_recall(s,block_seq(s,:)==3)));
    rem_typeblock0(s,4) = sum(num_recall(s,block_seq(s,:)==4)==0)/sum(~isnan(num_recall(s,block_seq(s,:)==4)));
    rem_typeblock0(s,5) = sum(num_recall(s,block_seq(s,:)==5)==0)/sum(~isnan(num_recall(s,block_seq(s,:)==5)));
    rem_typeblock0(s,6) = sum(num_recall(s,block_seq(s,:)==6)==0)/sum(~isnan(num_recall(s,block_seq(s,:)==6)));
    
    % rem1
    rem_typeblock1(s,1) = sum(num_recall(s,block_seq(s,:)==1)==1)/sum(~isnan(num_recall(s,block_seq(s,:)==1)));
    rem_typeblock1(s,2) = sum(num_recall(s,block_seq(s,:)==2)==1)/sum(~isnan(num_recall(s,block_seq(s,:)==2)));
    rem_typeblock1(s,3) = sum(num_recall(s,block_seq(s,:)==3)==1)/sum(~isnan(num_recall(s,block_seq(s,:)==3)));
    rem_typeblock1(s,4) = sum(num_recall(s,block_seq(s,:)==4)==1)/sum(~isnan(num_recall(s,block_seq(s,:)==4)));
    rem_typeblock1(s,5) = sum(num_recall(s,block_seq(s,:)==5)==1)/sum(~isnan(num_recall(s,block_seq(s,:)==5)));
    rem_typeblock1(s,6) = sum(num_recall(s,block_seq(s,:)==6)==1)/sum(~isnan(num_recall(s,block_seq(s,:)==6)));
    
    % rem2
    rem_typeblock2(s,1) = sum(num_recall(s,block_seq(s,:)==1)==2)/sum(~isnan(num_recall(s,block_seq(s,:)==1)));
    rem_typeblock2(s,2) = sum(num_recall(s,block_seq(s,:)==2)==2)/sum(~isnan(num_recall(s,block_seq(s,:)==2)));
    rem_typeblock2(s,3) = sum(num_recall(s,block_seq(s,:)==3)==2)/sum(~isnan(num_recall(s,block_seq(s,:)==3)));
    rem_typeblock2(s,4) = sum(num_recall(s,block_seq(s,:)==4)==2)/sum(~isnan(num_recall(s,block_seq(s,:)==4)));
    rem_typeblock2(s,5) = sum(num_recall(s,block_seq(s,:)==5)==2)/sum(~isnan(num_recall(s,block_seq(s,:)==5)));
    rem_typeblock2(s,6) = sum(num_recall(s,block_seq(s,:)==6)==2)/sum(~isnan(num_recall(s,block_seq(s,:)==6)));

end 
rem_typeblock0=rem_typeblock0.*100;
rem_typeblock1=rem_typeblock1.*100;
rem_typeblock2=rem_typeblock2.*100;

figure(31);clf;
for fig=1
    err = std(rem_typeblock2);
    b=bar(mean(rem_typeblock2,1),0.8);hold on
    er = errorbar(mean(rem_typeblock2,1),err);
    scatter(ones(size(rem_typeblock2(:,1))).*(1+(rand(size(rem_typeblock2(:,1)))-0.5)/5),rem_typeblock2(:,1),'k','filled');
    scatter(ones(size(rem_typeblock2(:,2))).*(2+(rand(size(rem_typeblock2(:,2)))-0.5)/5),rem_typeblock2(:,2),'k','filled');
    scatter(ones(size(rem_typeblock2(:,3))).*(3+(rand(size(rem_typeblock2(:,3)))-0.5)/5),rem_typeblock2(:,3),'k','filled');
    scatter(ones(size(rem_typeblock2(:,4))).*(4+(rand(size(rem_typeblock2(:,4)))-0.5)/5),rem_typeblock2(:,4),'k','filled');
    scatter(ones(size(rem_typeblock2(:,5))).*(5+(rand(size(rem_typeblock2(:,5)))-0.5)/5),rem_typeblock2(:,5),'k','filled');
    scatter(ones(size(rem_typeblock2(:,6))).*(6+(rand(size(rem_typeblock2(:,6)))-0.5)/5),rem_typeblock2(:,6),'k','filled');hold off
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 3;
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
    b.CData([1,4],:) = [0,0,1;0,0,1];
    b.CData([2,5],:) = [1,0,0;1,0,0];
    b.CData([3,6],:) = [0,1,0;0,1,0];
    xticks([1 2 3 4 5 6])
    xticklabels({'F-P-O','P-F-O','O-F-P','F-O-P','P-O-F','O-P-F'})
    ylim([0,100])
    ylabel('Percentage of trials')
    title('Image Recall 2')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end

% Separate by block number
rem_blocknumber = zeros(size(num_recall,1),6);
rem_blocknumber(:,1) = sum(num_recall(:,1:12)==2,2)./sum(~isnan(num_recall(:,1:12)),2)*100;
rem_blocknumber(:,2) = sum(num_recall(:,13:24)==2,2)./sum(~isnan(num_recall(:,13:24)),2)*100;
rem_blocknumber(:,3) = sum(num_recall(:,25:36)==2,2)./sum(~isnan(num_recall(:,25:36)),2)*100;
rem_blocknumber(:,4) = sum(num_recall(:,36:48)==2,2)./sum(~isnan(num_recall(:,36:48)),2)*100;
rem_blocknumber(:,5) = sum(num_recall(:,49:60)==2,2)./sum(~isnan(num_recall(:,49:60)),2)*100;
rem_blocknumber(:,6) = sum(num_recall(:,61:72)==2,2)./sum(~isnan(num_recall(:,61:72)),2)*100;

figure(32);clf;
for fig=1
    err = std(rem_blocknumber);
    b=bar(mean(rem_blocknumber,1),0.8);hold on
    er = errorbar(mean(rem_blocknumber,1),err);
    scatter(ones(size(rem_blocknumber(:,1))).*(1+(rand(size(rem_blocknumber(:,1)))-0.5)/5),rem_blocknumber(:,1),'k','filled');
    scatter(ones(size(rem_blocknumber(:,2))).*(2+(rand(size(rem_blocknumber(:,2)))-0.5)/5),rem_blocknumber(:,2),'k','filled');
    scatter(ones(size(rem_blocknumber(:,3))).*(3+(rand(size(rem_blocknumber(:,3)))-0.5)/5),rem_blocknumber(:,3),'k','filled');
    scatter(ones(size(rem_blocknumber(:,4))).*(4+(rand(size(rem_blocknumber(:,4)))-0.5)/5),rem_blocknumber(:,4),'k','filled');
    scatter(ones(size(rem_blocknumber(:,5))).*(5+(rand(size(rem_blocknumber(:,5)))-0.5)/5),rem_blocknumber(:,5),'k','filled');
    scatter(ones(size(rem_blocknumber(:,6))).*(6+(rand(size(rem_blocknumber(:,6)))-0.5)/5),rem_blocknumber(:,6),'k','filled');hold off
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 3;
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
%     b.CData([1,4],:) = [0,0,1;0,0,1];
%     b.CData([2,5],:) = [1,0,0;1,0,0];
%     b.CData([3,6],:) = [0,1,0;0,1,0];
    xticks([1 2 3 4 5 6])
    xticklabels({'Block1','Block2','Block3','Block4','Block5','Block6'})
    ylim([0,100])
    ylabel('Percentage of trials')
    title('Image Recall 2')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end

% Separate by item type
rem_item=zeros(size(num_recall,1),3);
for s = 1:size(num_recall,1)
    ind_fa = item_recall(s,enc_cat_index(s,:)==1);
    rem_item(s,1)=sum(ind_fa==1)/sum(ind_fa==1 |ind_fa==0)*100;
    
    ind_pl = item_recall(s,enc_cat_index(s,:)==2);
    rem_item(s,2)=sum(ind_pl==1)/sum(ind_pl==1 |ind_pl==0)*100;
    
    ind_ob = item_recall(s,enc_cat_index(s,:)==3);
    rem_item(s,3)=sum(ind_ob==1)/sum(ind_ob==1 |ind_ob==0)*100;
    clearvars ind_fa ind_pl ind_ob
end 


figure(4);clf;
for fig=1
    err = std(rem_item);
    b=bar(mean(rem_item,1),0.8);hold on
    er = errorbar(mean(rem_item,1),err);
    scatter(ones(size(rem_item(:,1))).*(1+(rand(size(rem_item(:,1)))-0.5)/5),rem_item(:,1),'k','filled');
    scatter(ones(size(rem_item(:,2))).*(2+(rand(size(rem_item(:,2)))-0.5)/5),rem_item(:,2),'k','filled');
    scatter(ones(size(rem_item(:,3))).*(3+(rand(size(rem_item(:,3)))-0.5)/5),rem_item(:,3),'k','filled');hold off
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 3;
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
    b.CData(1,:) = [0,0,1];
    b.CData(2,:) = [1,0,0];
    b.CData(3,:) = [0,1,0];
    xticks([1 2 3 4 5 6])
    xticklabels({'Face','Place','Object'})
    ylim([0,100])
    ylabel('Percentage of trials')
    title('Item recall')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end



% Median split
rem_num_recall = sum(num_recall==2,2)./sum(~isnan(num_recall),2)*100;
for_num_recall = sum(num_recall==0|num_recall==1,2)./sum(~isnan(num_recall),2)*100;

mean(rem_num_recall)
mean(for_num_recall)
std(rem_num_recall)
std(for_num_recall)

% wilcoxon
[p,h,stats] = signrank(rem_num_recall,for_num_recall);

num_ms=[rem_num_recall,for_num_recall];
figure(4);clf;
for fig=1
err = std(num_ms);
b=bar(mean(num_ms,1),0.8);hold on
er = errorbar(mean(num_ms,1),err); 
scatter(ones(size(num_ms(:,1))).*(1+(rand(size(num_ms(:,1)))-0.5)/5),num_ms(:,1),'k','filled');
scatter(ones(size(num_ms(:,2))).*(2+(rand(size(num_ms(:,2)))-0.5)/5),num_ms(:,2),'k','filled');hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;
b.FaceColor = 'flat';
b.FaceAlpha = 0.6;
    b.CData(1,:) = [0,0,1];
    b.CData(2,:) = [1,0,0];
xticks([1 2])
xticklabels({'Successful Recall','Unsuccessful Recall'})
ylim([0,100])
ylabel('Percentage of trials')
title('Image Recall')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end 


%% Difficulty rating
mean_diff = mean(diff_rating,2,'omitnan');
SD_diff= std(mean_diff);


% Percentage
mean_per_diff = zeros(30,4);
mean_per_diff(:,1)=sum(diff_rating==1,2)./sum(~isnan(diff_rating),2)*100;
mean_per_diff(:,2)=sum(diff_rating==2,2)./sum(~isnan(diff_rating),2)*100;
mean_per_diff(:,3)=sum(diff_rating==3,2)./sum(~isnan(diff_rating),2)*100;
mean_per_diff(:,4)=sum(diff_rating==4,2)./sum(~isnan(diff_rating),2)*100;

SD_per_diff=std(mean_per_diff);


[table,chi2,p] = crosstab(mean(num_recall,2,'omitnan'),mean(diff_rating,2,'omitnan'));

corrcoef(mean(num_recall,2,'omitnan'),mean(diff_rating,2,'omitnan'))


figure(5);clf;
for fig=1
err = std(mean_per_diff);
b=bar(mean(mean_per_diff,1),0.8);hold on
er = errorbar(mean(mean_per_diff,1),err); 
scatter(ones(size(mean_per_diff(:,1))).*(1+(rand(size(mean_per_diff(:,1)))-0.5)/5),mean_per_diff(:,1),'k','filled');
scatter(ones(size(mean_per_diff(:,2))).*(2+(rand(size(mean_per_diff(:,2)))-0.5)/5),mean_per_diff(:,2),'k','filled');
scatter(ones(size(mean_per_diff(:,3))).*(3+(rand(size(mean_per_diff(:,1)))-0.5)/5),mean_per_diff(:,1),'k','filled');
scatter(ones(size(mean_per_diff(:,4))).*(4+(rand(size(mean_per_diff(:,2)))-0.5)/5),mean_per_diff(:,2),'k','filled');hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;
b.FaceColor = 'flat';
%b.FaceAlpha = 0.6;
    %b.CData(1,:) = [0,0,1];
    %b.CData(2,:) = [1,0,0];
xticks([1 2 3 4])
xticklabels({'1','2','3','4'})
ylim([0,90])
ylabel('Percentage of trials')
xlabel('Difficulty Rating')
%title('Difficulty Rating')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end 

% difficulty rating after median-split
for s=1:30
    rem_diff(s) =  sum(diff_rating(s,num_recall(s,:)==2),2,'omitnan')./sum(~isnan(diff_rating(s,num_recall(s,:)==2)));
    for_diff(s) =  sum(diff_rating(s,num_recall(s,:)==1 | num_recall(s,:)==0),2,'omitnan')./sum(~isnan(diff_rating(s,num_recall(s,:)==1 | num_recall(s,:)==0)));
end 

mean(rem_diff)
mean(for_diff)
std(rem_diff)
std(for_diff)

[h,p,ci,stats] = ttest(rem_diff,for_diff);


diff_ms = [rem_diff',for_diff'];
indv_cond1 = ones(size(diff_ms(:,1))).*(1+(rand(size(diff_ms(:,1)))-0.5)/5);
indv_cond2 = ones(size(diff_ms(:,2))).*(2+(rand(size(diff_ms(:,2)))-0.5)/5);
C=repmat([1,0,0;0,0,1],1,1);
figure(6); clear clf;
boxplot(diff_ms,'Widths',0.6,'OutlierSize',4,'Symbol','');hold on
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),C(j,:),'FaceAlpha',.6);
end
set(findobj(gca,'type','line'),'linew',1.5)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
scatter(indv_cond1,diff_ms(:,1),'k','filled');
scatter(indv_cond2,diff_ms(:,2),'k','filled');
line([indv_cond1,indv_cond2]',[diff_ms(:,1),diff_ms(:,2)]','color',[0.5,0.5,0.5]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
% add sig line
line([1 2]',[3.2 3.2]','color',[0.5,0.5,0.5],'linewidth',3);
text([1.5], [3.3], 'n.s.','color',[0.5,0.5,0.5],'FontSize',14);hold off
ylabel('Difficulty rating')
xticks([1 2])
xticklabels({'Successful Recall','Unsuccessful Recall'})
ylim([1,3.5])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);