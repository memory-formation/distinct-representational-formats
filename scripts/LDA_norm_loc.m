%%
savepath = '';

for loaddata = 1
low_pass=0.5;  
savename = strcat(savepath,'ERP_loc_0.5Hzlowpass_group_Gaussian_smo50.mat');% 25 smo
load (savename);

savename = strcat(savepath,'Beh&index_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);

savename = strcat(savepath,'reject_ind_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);
end 
%% Normalization
for s = 1:length(ERP_loc_ds)
    disp(s)
    ind_loc = reject_index_loc(s,:)==1 & loc_hit_index(s,:)==1;
    
    data_loc = ERP_loc_ds(s).loc;

    ERP_loc_norm(s).data=zeros(size(data_loc));
    
    for t= 1:size(data_loc,2)
        for chan=1:62
            data2norm = squeeze(data_loc(chan,t,ind_loc));
            ERP_loc_norm(s).data(chan,t,ind_loc)=zscore(data2norm);
            clearvars data2norm
        end
    end
end



savename = strcat(savepath,'ERP_loc_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
save (savename,'ERP_loc_norm','-v7.3');


%% Start

for loaddata = 1

low_pass=0.5;  
savename = strcat(savepath,'ERP_loc_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
load (savename);

savename = strcat(savepath,'Beh&index_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);

savename = strcat(savepath,'reject_ind_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);
end 


tw = 1; % round 50 ms
for s=1:length(ERP_loc_norm)
    alldata =ERP_loc_norm(s).data;
    alldata([2,31,34],:,:)=[];% take out eyes channels.
    Xalldata = alldata(:,:,loc_hit_index(s,:) == 1 & reject_index_loc(s,:) ==1 );
    Yalltrain = loc_cat_index(s,loc_hit_index(s,:) == 1 & reject_index_loc(s,:) ==1);% labels
    
    for itime=1:size(Xalldata,2)-tw+1
        disp(['subject ' int2str(s) ' time '  int2str(itime) ])
        pred = zeros(1,size(Xalldata,3));
        dist = zeros(3,size(Xalldata,3));
        for itrial=1:size(Xalldata,3)
            t_feature = itime:itime+tw-1;
            Xtest = reshape(Xalldata(:,t_feature,itrial),1,[]); % just one trial to be classified
            Xtrain = Xalldata(:,t_feature,1:size(Xalldata,3)~=itrial);
            Xtrain = reshape(Xtrain,[size(Xtrain,1)*size(Xtrain,2),size(Xtrain,3)])';
            Ytrain = Yalltrain(1:size(Xalldata,3)~=itrial);
            
            [pred(itrial),~,dist(:,itrial),~] = xb_Multiclass_LDA(Xtest,Xtrain, Ytrain);
        end
        % general accu
        accuracy(s,itime)=sum(pred==Yalltrain)/size(Yalltrain,2);
        % each cat accuracy
        ind_fa = Yalltrain ==1;ind_pl = Yalltrain ==2;ind_ob = Yalltrain ==3;
        acc_fa(s,itime) = sum(pred(ind_fa)==Yalltrain(ind_fa))/sum(ind_fa);
        acc_pl(s,itime) = sum(pred(ind_pl)==Yalltrain(ind_pl))/sum(ind_pl);
        acc_ob(s,itime) = sum(pred(ind_ob)==Yalltrain(ind_ob))/sum(ind_ob);
        
        distance(s,itime).dis = sum(abs(dist),1);
    end
    
end

low_pass = 0.5;
tw = 1; % round 50 ms


savename = strcat(savepath,'LDA_loc_window(dist)_(Gaussian_norm)_',num2str(low_pass),'Hzlowpass_group_tw',int2str(tw), '_smo50.mat');% 25 smo
save (savename,'accuracy','distance','-v7.3');
savename = strcat(savepath,'LDA_loc_window(dist)_(Gaussian_norm)_eachcat_',num2str(low_pass),'Hzlowpass_group_tw',int2str(tw), '_smo50.mat');% 25 smo
save (savename,'acc_fa','acc_pl','acc_ob','-v7.3');
clearvars savepath savename


%% Plotting real
for plotting =1
plot_areaerrorbar(accuracy)
line([0,size(accuracy,2)],[0.33,0.33],'color','r','LineWidth', 2)
xticks([0 0.5 1 1.5 2 2.5]*512/5)
xticklabels({'0','500','1000','1500','2000','2500'})
xlabel('Encoding [ms]')
ylabel('Accuracy')
legend({'ACC','','Chance level'})
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);


dist_ave = zeros(size(distance));
for s = 1:size(distance,1)
    disp(int2str(s))
    for itime = 1:size(distance,2)
     dist_ave(s,itime)=mean(sum(abs(distance(s,itime).dis),1));
    end
end 

options.handle     = figure(30);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.alpha      = 0.7;
options.line_width = 2;
options.error      = 'sem';

plot_areaerrorbar(dist_ave,options)
xticks([0 0.5 1 1.5 2 2.5]*512/5)
xticklabels({'0','500','1000','1500','2000','2500'})
xlabel('Encoding [ms]')
ylabel('Distance')
title('Distance')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

end 

%% Permutation for surrogae trials 


for loaddata = 1

low_pass=0.5;  
savename = strcat(savepath,'ERP_loc_0.5Hzlowpass_group_Gaussian_smo50_withnorm.mat');% 25 smo
load (savename);

savename = strcat(savepath,'Beh&index_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);

savename = strcat(savepath,'reject_ind_loc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);
end 


tw = 1; % round 50 ms
n_perm = 50;
for i_perm=1:n_perm
    for s=1:length(ERP_loc_ds)
        disp(['subject ' int2str(s)  '/ perm ' int2str(i_perm)])
        alldata =ERP_loc_ds(s).loc;
        alldata([2,31,34],:,:)=[];% take out eyes channels.
        Xalldata = alldata(:,:,loc_hit_index(s,:) == 1 & reject_index_loc(s,:) ==1 );
        Yalltrain = loc_cat_index(s,loc_hit_index(s,:) == 1 & reject_index_loc(s,:) ==1);% labels
        
        Yalltrain_sh = Yalltrain(randperm(length(Yalltrain)));
        
        for itime=1:size(Xalldata,2)-tw+1
            
            pred = zeros(1,size(Xalldata,3));
            dist = zeros(3,size(Xalldata,3));
            for itrial=1:size(Xalldata,3)
                t_feature = itime:itime+tw-1;
                Xtest = reshape(Xalldata(:,t_feature,itrial),1,[]); % just one trial to be classified
                Xtrain = Xalldata(:,t_feature,1:size(Xalldata,3)~=itrial);
                Xtrain = reshape(Xtrain,[size(Xtrain,1)*size(Xtrain,2),size(Xtrain,3)])';
                Ytrain = Yalltrain_sh(1:size(Xalldata,3)~=itrial);
                
                [pred(itrial),~,dist(:,itrial),~] = xb_Multiclass_LDA(Xtest,Xtrain, Ytrain);
            end
            % general accu
            accuracy(s,itime,i_perm)=sum(pred==Yalltrain_sh)/size(Yalltrain_sh,2);
            % each cat accuracy
            ind_fa = Yalltrain_sh ==1;ind_pl = Yalltrain_sh ==2;ind_ob = Yalltrain_sh ==3;
            acc_fa(s,itime,i_perm) = sum(pred(ind_fa)==Yalltrain_sh(ind_fa))/sum(ind_fa);
            acc_pl(s,itime,i_perm) = sum(pred(ind_pl)==Yalltrain_sh(ind_pl))/sum(ind_pl);
            acc_ob(s,itime,i_perm) = sum(pred(ind_ob)==Yalltrain_sh(ind_ob))/sum(ind_ob);
            
            distance(s,itime,i_perm).dis = sum(abs(dist),1);
        end
        
    end
end
savepath = 'HERE PLACE YOUR PATH\'; 
savename = strcat(savepath,'LDA_loc_window(dist)_(Gaussian)_permpart1_',num2str(low_pass),'Hzlowpass_group_tw',int2str(tw), '.mat');% 25 smo
save (savename,'accuracy','distance','acc_fa','acc_pl','acc_ob','-v7.3');
