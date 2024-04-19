%% Normalize encoding data

savepath = '';
for loaddata = 1
    low_pass = 0.5; N = 50;

    
    savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
    load (savename);%save rejection index
    
    savename = strcat(savepath,'Beh&index_loc_',num2str(low_pass),'Hzlowpass_group.mat');
    load (savename);
    
    
    savename = strcat(savepath,'ERP_ds_(Gaussian)_',num2str(low_pass),'Hzlowpass_group_','smo_',num2str(N),'.mat');
    load (savename);%load ERPs
end

for s = 1:length(ERP_erp_ds)
    disp(s)
    ind_enc = reject_index_img(s,:)==1 ;
    
    data_enc = ERP_erp_ds(s).img;

    ERP_enc_norm(s).data=zeros(size(data_enc));
    
    for t= 1:size(data_enc,2)
        for chan=1:62
            data2norm = squeeze(data_enc(chan,t,ind_enc));
            ERP_enc_norm(s).data(chan,t,ind_enc)=zscore(data2norm);
            clearvars data2norm
        end
    end
end


savename = strcat(savepath,'ERP_enc_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
save (savename,'ERP_enc_norm','-v7.3');


%% Load normalized data
savepath = '';
for loaddata = 1

low_pass = 0.5; 
savename = strcat(savepath,'ERP_loc_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
load (savename);

savename = strcat(savepath,'Beh&index_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);

savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs

savename = strcat(savepath,'reject_ind_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);

savename = strcat(savepath,'ERP_enc_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
load (savename);%load ERPs
end

%% LDA
n_train = 13:23;
n_test  = 1:205;% all points

LDA_enc_prob_fa_ave = zeros(length(ERP_enc_norm),length(n_train),length(n_test),size(ERP_enc_norm(1).data,3));
LDA_enc_prob_pl_ave = zeros(length(ERP_enc_norm),length(n_train),length(n_test),size(ERP_enc_norm(1).data,3));
LDA_enc_prob_ob_ave = zeros(length(ERP_enc_norm),length(n_train),length(n_test),size(ERP_enc_norm(1).data,3));
LDA_enc_pred_ave    = zeros(length(ERP_enc_norm),length(n_train),length(n_test),size(ERP_enc_norm(1).data,3));

for t = n_train
    for s=1:length(ERP_enc_norm)
        disp(['subject_' int2str(s)  '  timewindow_' int2str(t)])
        testdata = ERP_enc_norm(s).data;
        testdata([2,31,34],:,:)=[];% take out eyes channels.
        
        locdata = ERP_loc_norm(s).data;
        locdata([2,31,34],:,:)=[];% take out eyes channels.
        
        traindata = mean(locdata(:,t, loc_hit_index(s,:) == 1 & reject_index_loc(s,:) ==1 ),2); % select the peak feature with time window from all clean and correct trials
        traindata = reshape(traindata,[size(traindata,1)*size(traindata,2),size(traindata,3)])';
        
        Ytrain = loc_cat_index(s,loc_hit_index(s,:) == 1 & reject_index_loc(s,:) ==1);% labels
        
        for itime=n_test
            
            i_testdata = testdata(:,itime,:);
            i_testdata = reshape(i_testdata,[size(i_testdata,1)*size(i_testdata,2),size(i_testdata,3)])';
            
            [pred,prob,dist,~] = xb_Multiclass_LDA(i_testdata,traindata, Ytrain);
            
            LDA_enc_prob_fa_ave(s,t-n_train(1)+1,itime,:) = prob(:,1);
            LDA_enc_prob_pl_ave(s,t-n_train(1)+1,itime,:) = prob(:,2);
            LDA_enc_prob_ob_ave(s,t-n_train(1)+1,itime,:) = prob(:,3);
            LDA_enc_pred_ave(s,t-n_train(1)+1,itime,:) = pred;
            
            distance{s,itime,t-n_train(1)+1} = sum(abs(dist));
            
        end
        
    end
end


savename = strcat(savepath,'LDA_predenc_p2p(+-50aroundpeak18)_withnorm', '.mat');
save (savename,'LDA_enc_prob_fa_ave','LDA_enc_prob_pl_ave','LDA_enc_prob_ob_ave','LDA_enc_pred_ave','distance','-v7.3');

%% Get averaged data
dist_S1 = zeros(size(distance));dist_S2 = zeros(size(distance));dist_S3 = zeros(size(distance));
for s = 1:size(distance,1)
    disp(s)
    ind_rej = reject_index_img(s,:)==1;
    ind_seq = rep_num(1:3,1,72); % always 1 2 3
    for t_loc=1:size(distance,3)
        for t_enc = 1:size(distance,2)
            dist = distance{s,t_enc,t_loc};
            dist_S1(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & ind_seq==1)),1));
            dist_S2(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & ind_seq==2)),1));
            dist_S3(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & ind_seq==3)),1));
            
            dist_all(s,t_enc,t_loc)= mean(sum(abs(dist(:,ind_rej )),1));

        end
    end

end 

savename = strcat(savepath,'LDA_predenc_avep2p(+-50aroundpeak18)_withnorm.mat');
save (savename,'dist_S1','dist_S2','dist_S3','dist_all','-v7.3');

% by category
dist_fa = zeros(size(distance));dist_pl = zeros(size(distance));dist_ob = zeros(size(distance));
for s = 1:size(distance,1)
    disp(s)
    ind_rej = reject_index_img(s,:)==1;
    ind_cat = enc_cat_index(s,:); % always 1 2 3
    for t_loc=1:size(distance,3)
        for t_enc = 1:size(distance,2)
            dist = distance{s,t_enc,t_loc};
            dist_fa(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & ind_cat==1)),1));
            dist_pl(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & ind_cat==2)),1));
            dist_ob(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & ind_cat==3)),1));
            
            %dist_all(s,t_enc,t_loc)= mean(sum(abs(dist(:,ind_rej )),1));

        end
    end

end 

%% Ploting 
figure(11);
plot(mean(mean(dist_all,3),1))


figure(122);
plot(mean(mean(dist_S1,3),1),'linewidth',2);hold on
plot(mean(mean(dist_S2,3),1),'linewidth',2);
plot(mean(mean(dist_S3,3),1),'linewidth',2);hold off
ylim([5.3,6.3])
xticks([0 0.5 1 1.5 2]*512/5)
xticklabels({'0','500','1000','1500','2000'})
xlabel('Encoding [ms]')
ylabel('Distance')
title('S1 S2 S3')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

% category exploration
for fig=1
figure(222);
plot(mean(mean(dist_fa,3),1),'linewidth',2);hold on
plot(mean(mean(dist_pl,3),1),'linewidth',2);
plot(mean(mean(dist_ob,3),1),'linewidth',2);hold off
ylim([5.3,6.3])
xticks([0 0.5 1 1.5 2]*512/5)
xticklabels({'0','500','1000','1500','2000'})
xlabel('Encoding [ms]')
ylabel('Distance')
title('Face Place Object')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

end 


for s=1:30
    ind_rej = reject_index_img(s,:)==1;
    ind_seq = rep_num(1:3,1,72); % always 1 2 3
    for t=1:11
        accu_S11(s,:,t) = sum(squeeze(LDA_enc_pred_ave(s,t,:,ind_rej & ind_seq==1)) == enc_cat_index(s,ind_rej & ind_seq==1),2)/sum(ind_rej & ind_seq==1);
        accu_S22(s,:,t) = sum(squeeze(LDA_enc_pred_ave(s,t,:,ind_rej & ind_seq==2)) == enc_cat_index(s,ind_rej & ind_seq==2),2)/sum(ind_rej & ind_seq==2);
        accu_S33(s,:,t) = sum(squeeze(LDA_enc_pred_ave(s,t,:,ind_rej & ind_seq==3)) == enc_cat_index(s,ind_rej & ind_seq==3),2)/sum(ind_rej & ind_seq==3);
    end
end


figure(13);
plot(mean(mean(accu_S11,3),1),'linewidth',2);hold on
plot(mean(mean(accu_S22,3),1),'linewidth',2);
plot(mean(mean(accu_S33,3),1),'linewidth',2);hold off
ylim([0.28,0.46])
xticks([0 0.5 1 1.5 2]*512/5)
xticklabels({'0','500','1000','1500','2000'})
xlabel('Encoding [ms]')
ylabel('Acc')
title('S1 S2 S3')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

% figure
for fig=13
    options.handle     = figure(fig);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(mean(dist_S1,3),options);hold on


options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;

plot_areaerrorbar(mean(dist_S2,3),options);hold on

options.color_area = [128 169 114]./255;    % Orange theme
options.color_line = [52  112  22]./255;

plot_areaerrorbar(mean(dist_S3,3),options);hold off

xticks([0 0.5 1 1.5 2 2.5]*512/5)
xticklabels({'0','500','1000','1500','2000','2500'})
xlabel('Encoding [ms]')
ylabel('Distance')
title('S1 S2 S3')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
end

%% Regression

for loaddata = 1
    low_pass = 0.5; N = 50;

    savename = strcat(savepath,'LDA_predenc_p2p(+-50aroundpeak18)_withnorm', '.mat');
    load (savename);
    
    savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
    load (savename);%load ERPs
    
    savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
    load (savename);%save rejection index
    
    clearvars loaddata
end



savename = strcat(savepath,'regression_predenc_p2p_notnorm_data_smo(Gaussian)_50.mat');
load (savename)
clearvars Ydata_alltp
for s=1:30
    disp(s)
    for i_t = 1:size(distance,3)
        data_s = cat(1,distance{s,:,i_t});
        
        ind_trial = reject_index_img(s,:)==1;
        data(i_t,:,:) = data_s(:,ind_trial);
         
    end
    data_all{s} = squeeze(mean(data,1));
    clearvars data
    
end
Ydata_alltp = cat(2,data_all{:});

r_data = zeros(0);
Ydata_alltp_ds=zeros(0);
N=20;
for i_trial = 1:size(Ydata_alltp,2)
    disp([int2str(i_trial) '/6247'])
    data = Ydata_alltp(:,i_trial);
    data = smoothdata(Ydata_alltp(:,i_trial),'movmean',N);
    r_data = downsample(data,1);
    Ydata_alltp_ds(:,i_trial)=r_data;
end

t_map_real = zeros(4,size(Ydata_alltp_ds,1)); %1,intercept; 2, Seq   3, Remfor 4. interaction
beta_map_real = zeros(4,size(Ydata_alltp_ds,1)); %1,intercept; 2, Remfor 
p_map_real = zeros(4,size(Ydata_alltp_ds,1)); %1,intercept; 2, Remfor 

% vif_result=zeros(2,2704);
tic
for i_time =1:size(Ydata_alltp_ds,1)
    if floor(i_time/20)==ceil(i_time/20)
        disp([int2str(i_time/size(Ydata_alltp_ds,1)*100) '% complete'])
    end
    
    tbl_time = table(double(Ydata_alltp_ds(i_time,:))',img_seq(:),subject_all(:),remfor(:)==1,'VariableNames',{'Dis','Img_seq','Subject_ID','RemFor'});
    lme_time = fitlme(tbl_time,'Dis~Img_seq+RemFor+Img_seq*RemFor+(1|Subject_ID)', 'CovariancePattern','Diagonal');
    
    
    
%     matrix_result = table2array(tbl_time);
%     vif_result(:,i_time) = vif(matrix_result(:,2:5));
    
    t_map_real(1,i_time) = lme_time.Coefficients{1,4};
    t_map_real(2,i_time) = lme_time.Coefficients{2,4};
    t_map_real(3,i_time) = lme_time.Coefficients{3,4};
    t_map_real(4,i_time) = lme_time.Coefficients{4,4};


    beta_map_real(1,i_time) = lme_time.Coefficients{1,2};
    beta_map_real(2,i_time) = lme_time.Coefficients{2,2};
    beta_map_real(3,i_time) = lme_time.Coefficients{3,2};
    beta_map_real(4,i_time) = lme_time.Coefficients{4,2};
    
    
    p_map_real(1,i_time) = lme_time.Coefficients{1,6};
    p_map_real(2,i_time) = lme_time.Coefficients{2,6};
    p_map_real(3,i_time) = lme_time.Coefficients{3,6};
    p_map_real(4,i_time) = lme_time.Coefficients{4,6};
    
    clearvars tbl_time lme_time
end
toc


tw_enc = 1:205; % first 2s data only

figure(3);
imagesc(t_map_real(2,tw_enc));
caxis([-2.5,2.5])

figure(4);
imagesc(p_map_real(2,tw_enc)<0.05/length(p_map_real(2,tw_enc)));

ind_sig = find(p_map_real(2,tw_enc)<0.05/length(p_map_real(2,tw_enc)));

% ploting
for fig=5
    for s=1:30
        LDA_s(:,s)=mean(Ydata_alltp_ds(tw_enc,subject_all(:)==s),2);
        
    end 
        
    
    
options.handle     = figure(fig);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(LDA_s',options);
xticks([0 0.5 1 1.5 2 ]*512/5)
xticklabels({'0','500','1000','1500','2000'})
xlabel('Encoding [ms]')
ylabel('Distance')
title('Distance')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
end

% ploting
for fig=6
    a=find(p_map_real(2,tw_enc)<0.05/length(p_map_real(2,tw_enc)));
    
figure(fig);clf;
plot(t_map_real(2,tw_enc),'linewidth',5);hold on
plot(t_map_real(3,tw_enc),'linewidth',5)
plot(t_map_real(4,tw_enc),'linewidth',5)
plot(a,repmat(3,length(a)),'.','MarkerSize',15,'color',[0.5,0.5,0.5]);
hold off
ylim([-5.5,5.5])
xticks([0 0.5 1 1.5 2 2.5 3 3.5]*512/5)
xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
xlabel('Encoding [ms]')
ylabel('t value')
title('Distance LMM')
legend({'Seq','Mem','Interaction'})
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
end


% control for accuracy
S1_acc_ave = mean(mean(accu_S11(:,ind_sig,:),3),2);
S2_acc_ave = mean(mean(accu_S22(:,ind_sig,:),3),2);
S3_acc_ave = mean(mean(accu_S33(:,ind_sig,:),3),2);

ave_acc_all = [S1_acc_ave,S2_acc_ave,S3_acc_ave];
figure(2);clf;
for fig=1
    err = std(ave_acc_all);
    b=bar(mean(ave_acc_all,1),0.8);hold on
    er = errorbar(mean(ave_acc_all,1),err);
    scatter(ones(size(ave_acc_all(:,1))).*(1+(rand(size(ave_acc_all(:,1)))-0.5)/5),ave_acc_all(:,1),'k','filled');
    scatter(ones(size(ave_acc_all(:,2))).*(2+(rand(size(ave_acc_all(:,2)))-0.5)/5),ave_acc_all(:,2),'k','filled');
    scatter(ones(size(ave_acc_all(:,3))).*(3+(rand(size(ave_acc_all(:,3)))-0.5)/5),ave_acc_all(:,3),'k','filled');
    %line([1 3]',[0.5 0.5]','color',[0.5,0.5,0.5],'linewidth',3);
    %text([2], [0.52], 'n.s.','color',[0.5,0.5,0.5],'FontSize',14);
    line([0.1 3.5]',[0.33 0.33]','LineStyle','--','color',[0.5,0.5,0.5],'linewidth',3);hold off
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 3;
    b.FaceColor = 'flat';
    b.FaceAlpha = 0.6;
    %b.CData([1,4],:) = [0,0,1;0,0,1];
    %b.CData([2,5],:) = [1,0,0;1,0,0];
    %b.CData([3,6],:) = [0,1,0;0,1,0];
    xticks([1 2 3 4 5 6])
    xticklabels({'S1','S2','S3'})
    ylim([0.25,0.5])
    ylabel('Numbers of items')
    title('Image Recall')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end

