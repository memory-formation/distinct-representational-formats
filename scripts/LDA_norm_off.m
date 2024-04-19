%% Normalize offset data
savepath = '';
for loaddata = 1
    low_pass = 0.5; N = 50;
    
    savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
    load (savename);%save rejection index
    
    savename = strcat(savepath,'ERP_enc_',num2str(low_pass),'Hzlowpass_group.mat');
    load (savename);
end

%% prepare data
% for prediction, we smooth and ds offset data
N= 50; % need to confirm the parameter
for subj= 1:size(ERP,2)
    disp(int2str(subj))
    data = ERP(subj).off;
    for ichan=1:62
        for itrial=1:size(data,3)
            data_smo(ichan,:,itrial)=smoothdata(data(ichan,:,itrial),'gaussian',N);
        end

    end
    for i=1:size(data_smo,3)
        %%%%%%%% enc %%%%%%%%%%%
        for j=1:62
            data_smo_down(j,:,i)=downsample(squeeze(data_smo(j,:,i)),5); % downsample and removal of the baseline
        end
    end
    clearvars i j
    ERP_off_ds(subj).data = data_smo_down;
    clearvars data_smo_down data_smo data
end 

savename = strcat(savepath,'ERP_ds_offset_(Gaussian)_',num2str(low_pass),'Hzlowpass_group_','smo_',num2str(N),'.mat');
save (savename,'ERP_off_ds','-v7.3');%load ERPs

%% Normalize data
for s = 1:length(ERP_off_ds)
    disp(s)
    ind_off = reject_index_off(s,:)==1 ;
    
    data_off = ERP_off_ds(s).data;

    ERP_off_norm(s).data=zeros(size(data_off));
    
    for t= 1:size(data_off,2)
        for chan=1:62
            data2norm = squeeze(data_off(chan,t,ind_off));
            ERP_off_norm(s).data(chan,t,ind_off)=zscore(data2norm);
            clearvars data2norm
        end
    end
end


savename = strcat(savepath,'ERP_off_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
save (savename,'ERP_off_norm','-v7.3');

%% LDA

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

savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save rejection index  

savename = strcat(savepath,'ERP_off_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
load (savename);%load ERPs
end

%% LDA
n_train = 1:205;
n_test  = 1:359;% all points

LDA_off_prob_fa_ave = zeros(length(ERP_off_norm),length(n_train),length(n_test),size(ERP_off_norm(1).data,3));
LDA_off_prob_pl_ave = zeros(length(ERP_off_norm),length(n_train),length(n_test),size(ERP_off_norm(1).data,3));
LDA_off_prob_ob_ave = zeros(length(ERP_off_norm),length(n_train),length(n_test),size(ERP_off_norm(1).data,3));
LDA_off_pred_ave    = zeros(length(ERP_off_norm),length(n_train),length(n_test),size(ERP_off_norm(1).data,3));

for t = n_train
    for s=1:length(ERP_off_norm)
        disp(['subject_' int2str(s)  '  timewindow_' int2str(t)])
        testdata = ERP_off_norm(s).data;
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
            
            LDA_off_prob_fa_ave(s,t-n_train(1)+1,itime,:) = prob(:,1);
            LDA_off_prob_pl_ave(s,t-n_train(1)+1,itime,:) = prob(:,2);
            LDA_off_prob_ob_ave(s,t-n_train(1)+1,itime,:) = prob(:,3);
            LDA_off_pred_ave(s,t-n_train(1)+1,itime,:) = pred;
            
            distance{s,itime,t-n_train(1)+1} = sum(abs(dist));
            
        end
        
    end
end



%savename = strcat(savepath,'LDA_predoff_p2p(+-50aroundpeak18)_withnorm', '.mat');
savename = strcat(savepath,'LDA_predoff_p2p(all)_withnorm', '.mat');
save (savename,'LDA_off_prob_fa_ave','LDA_off_prob_pl_ave','LDA_off_prob_ob_ave','LDA_off_pred_ave','distance','-v7.3');

%%
dist_high = zeros(size(distance));dist_low = zeros(size(distance));dist_all = zeros(size(distance));
for s = 1:size(distance,1)
    disp(s)
    ind_rej = reject_index_off(s,:)==1;

    for t_loc=1:size(distance,3)
        for t_enc = 1:size(distance,2)
            dist = distance{s,t_enc,t_loc};
            dist_high(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & num_recall(s,:)==2)),1));
            dist_low(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & (num_recall(s,:)==1 | num_recall(s,:)==0) )),1));
            
            dist_all(s,t_enc,t_loc)= mean(sum(abs(dist(:,ind_rej )),1));

        end
    end

end 

savename = strcat(savepath,'LDA_predoff_avep2p(+-50aroundpeak18)_withnorm.mat');
save (savename,'dist_high','dist_low','dist_all','-v7.3');

N=20;
for s = 1:size(distance,1)
    for i_t = 1:size(distance,3)
    
        data_s =  cat(1,distance{s,:,i_t});
        data_high(i_t,:,:) = data_s(:,reject_index_off(s,:)==1 & num_recall(s,:)==2);
        data_low(i_t,:,:) = data_s(:,reject_index_off(s,:)==1 &(num_recall(s,:)==1 | num_recall(s,:)==0));
    end
    data_avehigh = squeeze(mean(data_high,1));
    data_avelow  = squeeze(mean(data_low,1));
    
    for itrial = 1:size(data_avehigh,2)
        data_smohigh(:,itrial) = smoothdata(data_avehigh(:,itrial),'movmean',N);  
    end 
    for itrial = 1:size(data_avelow,2)
        data_smolow(:,itrial) = smoothdata(data_avelow(:,itrial),'movmean',N);  
    end 
    
    dist_all_high(s,:) = mean(data_smohigh,2);
    dist_all_low(s,:) = mean(data_smolow,2);
    clearvars data_*
end 


figure(1);clf;
plot(mean(dist_all_high,1),'linewidth',2);hold on
plot(mean(dist_all_low,1),'linewidth',2);hold off
legend({'High','Low'})
xticks([0 0.5 1 1.5 2 2.5 3 3.5]*512/5)
xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
xlabel('Offset [ms]')
ylabel('D')
ylim([5.5,6])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

[clusters, p_values, t_sums, permutation_distribution ] = permutest( dist_all_high',  dist_all_low', true, ...
    0.05, 1000, true );

[a,b,c,d] = ttest(mean(dist_high,3), mean(dist_low,3));


% ploting
for fig=30

    
options.handle     = figure(fig);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(dist_all_high,options);hold on


options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;

plot_areaerrorbar(dist_all_low,options);hold on
hold off

xticks([0 0.5 1 1.5 2 2.5 3 3.5]*512/5)
xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
xlabel('Encoding [ms]')
ylabel('Distance')
title('Distance')
legend({'Rem','','For'})
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
ylim([5,6.5])
end


%%  2D check

dist_high = zeros(size(distance));dist_low = zeros(size(distance));dist_all = zeros(size(distance));
for s = 1:size(distance,1)
    disp(s)
    ind_rej = reject_index_off(s,:)==1;

    for t_loc=1:size(distance,3)
        for t_enc = 1:size(distance,2)
            dist = distance{s,t_enc,t_loc};
            dist_high(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & num_recall(s,:)==2)),1));
            dist_low(s,t_enc,t_loc) = mean(sum(abs(dist(:,ind_rej & (num_recall(s,:)==1 | num_recall(s,:)==0) )),1));
            
            dist_all(s,t_enc,t_loc)= mean(sum(abs(dist(:,ind_rej )),1));

        end
    end

end 
figure(11);clf;
imagesc(squeeze(mean(dist_high,1))'); axis xy
figure(12);clf;
imagesc(squeeze(mean(dist_low,1))'); axis xy


[clusters, p_values, t_sums, permutation_distribution ] = permutest( permute(dist_high,[2,3,1]),  permute(dist_low,[2,3,1]), true, ...
    0.05, 1000, true );

