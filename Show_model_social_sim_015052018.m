function [seq] = Show_model_social_sim_015052018(params)
%SOCIAL_SIM Summary of this function goes here
%   Detailed explanation goes here

% model params:
%     [anxiety_0, anxiety_tau , reward_0_1, reward_0_2, reward_tau,
%     beta1_reward, beta2_reward, beta1_anxiety, beta2_anxiety,
%     choose_Stillness_score_0, choose_Exploration_score_0];

%%%% Mice SP model parameters(after optimization):
%%%% [7.66,10355.62,4.98,4.72,100406.23,1.10,1.12,-0.08,0.71,2.72,5.19]

%%%% Mice SNP model parameters options - same as mice SP but different reward values
%%%% [7.66,10355.62,4.7,4.9,100406.23,1.10,1.12,-0.08,0.71,2.72,5.19]
%%%% I am not sure it is needed since it is very similar to SP results

%%%% Rats SP model parameters options - same as mice SP but different reward values
%%%% [7.66,10355.62,5.8,4.72,33406.23,1.10,1.12,-0.08,0.71,2.72,5.19]

%%%% Rats SNP model parameters options - same as mice SP but different reward values
%%%% [7.66,10355.62,5.17,5.3,33406.23,1.10,1.12,-0.08,0.71,2.72,5.19]


   sim_length = 9000;  % Simulation duration
   num_anim = 60;      % Number of animals to simulate
   anxiety_0 = params(1); % Initial anxiety  
   anxiety_tau = params(2); % Anxiety decay time constant
   reward_0 = params(3:4); % Reward for each stimulus 
   reward_tau = params(5); % Reward time constant
   beta1_reward = params(6); % The reward effect on the decision to continue interaction
   beta2_reward = params(7);  % The reward effect on the decision to choose a particular stimulus
   beta1_anxiety = params(8); % The anxiety effect on the decision to continue interaction (should be negative)
   beta2_anxiety = params(9); % The anxiety effect on the decision to continue grooming
   choose_Stillness_score_0 = params(10);
   choose_Exploration_score_0 = params(11); 

   init_state = 4;         % Initial state. 1,2=stimuli, 3=choice position, 4=grooming

   seq = zeros(num_anim, sim_length);

   for anim = 1:num_anim
      s = init_state;
      anxiety = zeros(1, sim_length);
      reward = zeros(sim_length, 2);
      tot = zeros(1, 4);
      for t=1:sim_length
         seq(anim, t) = s;
         tot(s) = tot(s) + 1;
         
         %%%%%%%%%% Please change the parameters per animal:
         reward(t, 1) = reward_0(1)* exp(-tot(1)/reward_tau);      %%%% For rats                   
         reward(t, 2) = reward_0(2)* exp(-tot(2)/reward_tau);      %%%% For rats
%          reward(t, 1) = reward_0(1)* exp(tot(1)/reward_tau);      %%%% For mice                 
%          reward(t, 2) = reward_0(2)* exp(tot(2)/reward_tau);      %%%% For mice  
         
         anxiety(t) = anxiety_0 * exp(-t/anxiety_tau)-abs(reward(t, 1)-reward(t, 2));   %%%% Now depends also on Delta-rewards.            
         stay_stim_score = beta1_reward*reward(t,:) + beta1_anxiety*anxiety(t);
         choose_stim_score = beta2_reward*reward(t,:);
         choose_Stillness_score=choose_Stillness_score_0-beta2_anxiety*anxiety(t);
         choose_Exploration_score=choose_Exploration_score_0;
         stay_Exploration_score = choose_Stillness_score+choose_Exploration_score;
        
         stay_stim_prob = 1./(exp(-stay_stim_score)+1);
         choose_prob = softmax([choose_stim_score, choose_Stillness_score, choose_Exploration_score]')';
         stay_Exploration_prob = 1./(exp(-stay_Exploration_score)+1);
        
         p = zeros(4, 4);
         p(1, 1) = stay_stim_prob(1);
         p(1, 3) = 1 - stay_stim_prob(1);
         p(2, 2) = stay_stim_prob(2);
         p(2, 3) = 1 - stay_stim_prob(2);
         p(3, :) = choose_prob;
         p(4, 4) = stay_Exploration_prob;
         p(4, 3) = 1 - stay_Exploration_prob;
         s = new_state(p(s, :)); 
      end
   end

   %%%%%% Output data for model evaluation:
   [ShortLongInteractionsStimulus1,ShortLongInteractionsStimulus2, ... 
   TotalTimeStim1_Short_Interactions_AlongTime,TotalTimeStim2_Short_Interactions_AlongTime, ...
   TotalTimeStim1_Long_Interactions_AlongTime,TotalTimeStim2_Long_Interactions_AlongTime, ...
   TransitionsAlongTimeInMinutes_Minute1] = Evaluate(seq,num_anim,sim_length);

end

function [state] = new_state(probs)

   dice = rand(); 
   if dice<probs(1)
      state = 1;
   elseif dice>=probs(1) && dice<(probs(1)+probs(2))
      state = 2;
   elseif  dice>=(probs(1)+probs(2)) && dice<(probs(1)+probs(2)+probs(3))
      state = 3;
   elseif dice>=(probs(1)+probs(2)+probs(3))
      state = 4;
   end
   
end

function [ShortLongInteractionsStimulus1,ShortLongInteractionsStimulus2, ...
          TotalTimeStim1_Short_Interactions_AlongTime,TotalTimeStim2_Short_Interactions_AlongTime, ...
          TotalTimeStim1_Long_Interactions_AlongTime,TotalTimeStim2_Long_Interactions_AlongTime, ...
          TransitionsAlongTimeInMinutes_Minute1] = Evaluate(seq,num_anim,sim_length)

   Stimulus1TotalExplorationPopulationSummary=[];
   Stimulus2TotalExplorationPopulationSummary=[];
   Stimulus1ExplorationAlongTimePopulationSummary=[];
   Stimulus2ExplorationAlongTimePopulationSummary=[];
   Stimulus1DurationsOfInteractions=[];
   Stimulus2DurationsOfInteractions=[];
   for anim=1:num_anim
      %%%%%%%% Stimulus 1 %%%%%%%
      StimulusExplorationAlongTime=[];
      StimulusExplorationAlongTime(find(seq(anim,:)==1))=1;
      StimulusExplorationAlongTime=[0,StimulusExplorationAlongTime,zeros(1,sim_length+1-length(StimulusExplorationAlongTime))]; 
      InteractionStartEndTimes=[];
      InteractionStartEndTimes=diff([0,StimulusExplorationAlongTime,0]);
      
      %%%%% population summary of stimulus 1 exploration total time
      Stimulus1TotalExplorationPopulationSummary(1,anim)=length(find(StimulusExplorationAlongTime)); 
      %%%%% population summary of stimulus 1 exploration along time
      Stimulus1ExplorationAlongTimePopulationSummary=[Stimulus1ExplorationAlongTimePopulationSummary;StimulusExplorationAlongTime];
      %%%%% population summary of stimulus 1 number of interactions and duration of interactions entire time
      Stimulus1DurationsOfInteractions{1,anim}=find(InteractionStartEndTimes==-1)-find(InteractionStartEndTimes==1);
      
      %%%%%%%% Stimulus 2 %%%%%%%
      StimulusExplorationAlongTime=[];
      StimulusExplorationAlongTime(find(seq(anim,:)==2))=1;
      StimulusExplorationAlongTime=[0,StimulusExplorationAlongTime,zeros(1,sim_length+1-length(StimulusExplorationAlongTime))]; 
      InteractionStartEndTimes=[];
      InteractionStartEndTimes=diff([0,StimulusExplorationAlongTime,0]);
      
      %%%%% population summary of stimulus 2 exploration total time
      Stimulus2TotalExplorationPopulationSummary(1,anim)=length(find(StimulusExplorationAlongTime)); 
      %%%%% population summary of stimulus 2 exploration along time
      Stimulus2ExplorationAlongTimePopulationSummary=[Stimulus2ExplorationAlongTimePopulationSummary;StimulusExplorationAlongTime];
      %%%%% population summary of stimulus 2 number of interactions and duration of interactions entire time
      Stimulus2DurationsOfInteractions{1,anim}=find(InteractionStartEndTimes==-1)-find(InteractionStartEndTimes==1);
   end
   
   figure;

    Bins=300:600:length(Stimulus1ExplorationAlongTimePopulationSummary)-300;
    TransitionTimesForHistAll=[];
    TransitionTimesForHistAllForRaster={};
    Stimulus1OnlyIntervalsTimeHistAllAnimal=[];
    Stimulus2OnlyIntervalsTimeHistAllAnimal=[];
    Stimulus1EpoksTimesDurationsBinsAll=[];
    Stimulus2EpoksTimesDurationsBinsAll=[];
    Stimulus1IntervalTimesDurationsBinsAll=[];
    Stimulus2IntervalTimesDurationsBinsAll=[];
    Stimulus1InteractionHeattMap=Stimulus1ExplorationAlongTimePopulationSummary;
    Stimulus2InteractionHeattMap=Stimulus2ExplorationAlongTimePopulationSummary;
 
    for i=1:size(Stimulus1ExplorationAlongTimePopulationSummary,1)
   
       Stimulus1EpoksTimesDurationsBins(1:1:round(length(Stimulus1ExplorationAlongTimePopulationSummary)/1800),1:1:4)=0;
       Stimulus2EpoksTimesDurationsBins(1:1:round(length(Stimulus2ExplorationAlongTimePopulationSummary)/1800),1:1:4)=0; 
       Stimulus1Epoks=[];
       Stimulus1OnlyIntervalsSingleAnimal=[]; %%% Intervals after stimulus 1 which takes into account only epokes of stimulus 1 (the intervals include stimulus 2 epokes)
       Stimulus1InteractionsStart=[]; 
       Stimulus1InteractionsEnd=[];
       Stimulus2Epoks=[];
       Stimulus2OnlyIntervalsSingleAnimal=[];  %%% Intervals after stimulus 2 which takes into account only epokes of stimulus 2 (the intervals include stimulus 1 epokes)
       Stimulus1IntervalTimesDurationsBins(1:1:round(length(Stimulus1ExplorationAlongTimePopulationSummary)/1800),1:1:3)=0;
       Stimulus2IntervalTimesDurationsBins(1:1:round(length(Stimulus2ExplorationAlongTimePopulationSummary)/1800),1:1:3)=0; 
       Stimulus2InteractionsStart=[]; 
       Stimulus2InteractionsEnd=[];
   
       Stimulus1InteractionsStart=find(diff(Stimulus1ExplorationAlongTimePopulationSummary(i,:))==1);
       Stimulus1InteractionsEnd=find(diff(Stimulus1ExplorationAlongTimePopulationSummary(i,:))==-1);
       if ~isempty(Stimulus1InteractionsStart) && ~isempty(Stimulus1InteractionsEnd)
          if Stimulus1InteractionsStart(1)>Stimulus1InteractionsEnd(1)
             Stimulus1InteractionsEnd=Stimulus1InteractionsEnd(2:end);
          end
          if Stimulus1InteractionsStart(end)>Stimulus1InteractionsEnd(end)
             Stimulus1InteractionsStart=Stimulus1InteractionsStart(1:end-1);
          end
       end
       Stimulus2InteractionsStart=find(diff(Stimulus2ExplorationAlongTimePopulationSummary(i,:))==1);
       Stimulus2InteractionsEnd=find(diff(Stimulus2ExplorationAlongTimePopulationSummary(i,:))==-1);
       if ~isempty(Stimulus2InteractionsStart) && ~isempty(Stimulus2InteractionsEnd)
          if Stimulus2InteractionsStart(1)>Stimulus2InteractionsEnd(1)
             Stimulus2InteractionsEnd=Stimulus2InteractionsEnd(2:end);
          end
          if Stimulus2InteractionsStart(end)>Stimulus2InteractionsEnd(end)
             Stimulus2InteractionsStart=Stimulus2InteractionsStart(1:end-1);
          end
       end
   
       %%%%%%%%%% Interactions with different bouts durations along time %%%%%
   
       Stimulus1Epoks=(Stimulus1InteractionsStart(1:end)-Stimulus1InteractionsEnd(1:end))*-1;
       Stimulus2Epoks=(Stimulus2InteractionsStart(1:end)-Stimulus2InteractionsEnd(1:end))*-1;
       for Stimulus1EpoksNum=1:length(Stimulus1Epoks)
          TempStimulus1EpokeTimeHist=[];
          TempStimulus1EpokeDurationHist=[];
          TempStimulus1EpokeTimeHist=histcounts(Stimulus1InteractionsStart(Stimulus1EpoksNum),0:1800:size(Stimulus1ExplorationAlongTimePopulationSummary,2));
          TempStimulus1EpokeDurationHist=histcounts(Stimulus1Epoks(Stimulus1EpoksNum)/30,[0 6 19 120]);
          Stimulus1EpoksTimesDurationsBins(find(TempStimulus1EpokeTimeHist),find(TempStimulus1EpokeDurationHist))=Stimulus1EpoksTimesDurationsBins(find(TempStimulus1EpokeTimeHist),find(TempStimulus1EpokeDurationHist))+Stimulus1Epoks(Stimulus1EpoksNum)/30;    
       end
       for Stimulus2EpoksNum=1:length(Stimulus2Epoks)
          TempStimulus2EpokeTimeHist=[];
          TempStimulus2EpokeDurationHist=[];
          TempStimulus2EpokeTimeHist=histcounts(Stimulus2InteractionsStart(Stimulus2EpoksNum),0:1800:length(Stimulus1ExplorationAlongTimePopulationSummary));
          TempStimulus2EpokeDurationHist=histcounts(Stimulus2Epoks(Stimulus2EpoksNum)/30,[0 6 19 120]);
          Stimulus2EpoksTimesDurationsBins(find(TempStimulus2EpokeTimeHist),find(TempStimulus2EpokeDurationHist))=Stimulus2EpoksTimesDurationsBins(find(TempStimulus2EpokeTimeHist),find(TempStimulus2EpokeDurationHist))+Stimulus2Epoks(Stimulus2EpoksNum)/30;    
       end
   
       Stimulus1EpoksTimesDurationsBinsAll(:,:,i)=Stimulus1EpoksTimesDurationsBins;
       Stimulus2EpoksTimesDurationsBinsAll(:,:,i)=Stimulus2EpoksTimesDurationsBins;
 
       %%%%%%%%%%%%%%% Heat Map for bouts durations along time %%%%

       for Stimulus1EpoksNum=1:length(Stimulus1Epoks)
          Stimulus1InteractionHeattMap(i,Stimulus1InteractionsStart(Stimulus1EpoksNum):Stimulus1InteractionsEnd(Stimulus1EpoksNum))=Stimulus1Epoks(Stimulus1EpoksNum)/30;
       end
       for Stimulus2EpoksNum=1:length(Stimulus2Epoks)
          Stimulus2InteractionHeattMap(i,Stimulus2InteractionsStart(Stimulus2EpoksNum):Stimulus2InteractionsEnd(Stimulus2EpoksNum))=Stimulus2Epoks(Stimulus2EpoksNum)/30;
       end
    
       %%%%%%%%%%% Transitions %%%%%%%%%%%%%%%%%%%%%%%%
    
       LastStimulusToExplore=0;
       TransitionTimes=[];
       while min([Stimulus1InteractionsStart Stimulus2InteractionsStart])<100000
          StimulusToExplore=min([Stimulus1InteractionsStart Stimulus2InteractionsStart]);
          if find(Stimulus1InteractionsStart==StimulusToExplore)
             if LastStimulusToExplore==2
                TransitionTimes=[TransitionTimes,StimulusToExplore]; 
                LastStimulusToExplore=1;
                Stimulus1InteractionsStart(find(Stimulus1InteractionsStart==StimulusToExplore))=100001;
             elseif LastStimulusToExplore==0
                LastStimulusToExplore=1; 
                Stimulus1InteractionsStart(find(Stimulus1InteractionsStart==StimulusToExplore))=100001;
             else
                Stimulus1InteractionsStart(find(Stimulus1InteractionsStart==StimulusToExplore))=100001; 
             end
          elseif find(Stimulus2InteractionsStart==StimulusToExplore)
             if LastStimulusToExplore==1
                TransitionTimes=[TransitionTimes,StimulusToExplore]; 
                LastStimulusToExplore=2;
                Stimulus2InteractionsStart(find(Stimulus2InteractionsStart==StimulusToExplore))=100001;
             elseif LastStimulusToExplore==0
                LastStimulusToExplore=2; 
                Stimulus2InteractionsStart(find(Stimulus2InteractionsStart==StimulusToExplore))=100001;
             else
                Stimulus2InteractionsStart(find(Stimulus2InteractionsStart==StimulusToExplore))=100001; 
             end   
          end
       end
 
       TempTransitionTimesForHist=[];
       [TempTransitionTimesForHist,NotForUse] = hist(TransitionTimes,Bins);
       TransitionTimesForHistAll=[TransitionTimesForHistAll; TempTransitionTimesForHist];
       TransitionTimesForHistAllForRaster{i}=TransitionTimes;
    end

    ShortLongInteractionsStimulus1=[squeeze(sum(Stimulus1EpoksTimesDurationsBinsAll(:,1,:),1)) squeeze(sum(Stimulus1EpoksTimesDurationsBinsAll(:,2,:),1))  squeeze(sum(Stimulus1EpoksTimesDurationsBinsAll(:,3,:),1))  Stimulus1TotalExplorationPopulationSummary'/30];
    ShortLongInteractionsStimulus2=[squeeze(sum(Stimulus2EpoksTimesDurationsBinsAll(:,1,:),1)) squeeze(sum(Stimulus2EpoksTimesDurationsBinsAll(:,2,:),1))  squeeze(sum(Stimulus2EpoksTimesDurationsBinsAll(:,3,:),1))  Stimulus2TotalExplorationPopulationSummary'/30];
 
    TransitionsAlongTimeInMinutes=[];
    for i=1:size(TransitionTimesForHistAll,2)
       if mod(i,3)==1 && i+2<=size(TransitionTimesForHistAll,2)
          TransitionsAlongTimeInMinutes=[TransitionsAlongTimeInMinutes, mean(TransitionTimesForHistAll(:,i:i+2),2)];
       end
    end

    
    
    
    
    
    clear Stimulus1ExplorationAlongTimePopulationSummaryBined Stimulus2ExplorationAlongTimePopulationSummaryBined
    BinSize=600; %%% size of bin (num of frames) for exploration along session anaylsis.
    Stimulus1ExplorationAlongTimePopulationSummaryBined(1:1:size(Stimulus1ExplorationAlongTimePopulationSummary,1),1:1:floor(size(Stimulus1ExplorationAlongTimePopulationSummary,2)/600))=0;
    for l=1:size(Stimulus1ExplorationAlongTimePopulationSummary,1)
       for i=1:size(Stimulus1ExplorationAlongTimePopulationSummary,2)
          if mod(i,BinSize)==0
             Stimulus1ExplorationAlongTimePopulationSummaryBined(l,i/600)=sum(Stimulus1ExplorationAlongTimePopulationSummary(l,(i-BinSize+1):i),2);
          end
       end
    end
    Stimulus2ExplorationAlongTimePopulationSummaryBined(1:1:size(Stimulus2ExplorationAlongTimePopulationSummary,1),1:1:floor(size(Stimulus2ExplorationAlongTimePopulationSummary,2)/600))=0;
    for l=1:size(Stimulus2ExplorationAlongTimePopulationSummary,1)
       for i=1:size(Stimulus2ExplorationAlongTimePopulationSummary,2)
          if mod(i,BinSize)==0
             Stimulus2ExplorationAlongTimePopulationSummaryBined(l,i/600)=sum(Stimulus2ExplorationAlongTimePopulationSummary(l,(i-BinSize+1):i),2);
          end
       end
    end

    subplot(3,3,1); 
    hold on;
    errorbar(20:20:size(Stimulus1ExplorationAlongTimePopulationSummaryBined,2)*20,mean(Stimulus1ExplorationAlongTimePopulationSummaryBined/30,1),...
       std(Stimulus1ExplorationAlongTimePopulationSummaryBined/30,0,1)/sqrt(size(Stimulus1TotalExplorationPopulationSummary/30,2)),'g')
    errorbar(20:20:size(Stimulus1ExplorationAlongTimePopulationSummaryBined,2)*20,mean(Stimulus2ExplorationAlongTimePopulationSummaryBined/30,1),...
       std(Stimulus2ExplorationAlongTimePopulationSummaryBined/30,0,1)/sqrt(size(Stimulus2TotalExplorationPopulationSummary/30,2)),'r')
    xlim([0 size(Stimulus1ExplorationAlongTimePopulationSummaryBined,2)*20+10]);
    ylim([0 20]);
    xlabel('Time (20 Sec bins)')
    ylabel('Investigation time (sec)')
    title('Stimuli investigation along time (30Hz recording)')
    hold off;
    
    subplot(3,3,2); 
    hold on;
    bar([1 3 5 7],mean(ShortLongInteractionsStimulus1,1),0.4,'g');
    bar([2 4 6 8],mean(ShortLongInteractionsStimulus2,1),0.4,'r');     
    errorbar([1 3 5 7],mean(ShortLongInteractionsStimulus1,1),...
       std(ShortLongInteractionsStimulus1,0,1)/sqrt(size(ShortLongInteractionsStimulus1,1)),'g','LineStyle','none')
    errorbar([2 4 6 8],mean(ShortLongInteractionsStimulus2,1),...
       std(ShortLongInteractionsStimulus2,0,1)/sqrt(size(ShortLongInteractionsStimulus2,1)),'r','LineStyle','none')
    xlim([0 9]);
    ax = gca;
    ax.XTick=[1 3 5 7];
    ax.XTickLabel = ({'Short','Intermediate','Long','Total time'});
    ylabel('Time of investigations');
    xlabel('Bout duration')
    title('Time of different bouts');
    legend('Stimulus 1','Stimulus 2')
    hold off; 
    
    subplot(3,3,3); 
    hold on;
    bar([1:2:size(Stimulus1EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus1EpoksTimesDurationsBinsAll(:,1,:),3),0.4,'g');
    bar([2:2:size(Stimulus2EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus2EpoksTimesDurationsBinsAll(:,1,:),3),0.4,'r');
    errorbar([1:2:size(Stimulus1EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus1EpoksTimesDurationsBinsAll(:,1,:),3),std(Stimulus1EpoksTimesDurationsBinsAll(:,1,:),0,3)/sqrt(size(Stimulus1EpoksTimesDurationsBinsAll(:,1,:),3)),'g','LineStyle','none');
    errorbar([2:2:size(Stimulus2EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus2EpoksTimesDurationsBinsAll(:,1,:),3),std(Stimulus2EpoksTimesDurationsBinsAll(:,1,:),0,3)/sqrt(size(Stimulus2EpoksTimesDurationsBinsAll(:,1,:),3)),'r','LineStyle','none');
    xlabel('Time')
    ylabel('Investigation time')
    title('Short bouts along time')
    hold off; 

    subplot(3,3,4); 
    hold on;
    bar([1:2:size(Stimulus1EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus1EpoksTimesDurationsBinsAll(:,3,:),3),0.4,'g');
    bar([2:2:size(Stimulus2EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus2EpoksTimesDurationsBinsAll(:,3,:),3),0.4,'r');
    errorbar([1:2:size(Stimulus1EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus1EpoksTimesDurationsBinsAll(:,3,:),3),std(Stimulus1EpoksTimesDurationsBinsAll(:,3,:),0,3)/sqrt(size(Stimulus1EpoksTimesDurationsBinsAll(:,3,:),3)),'g','LineStyle','none');
    errorbar([2:2:size(Stimulus2EpoksTimesDurationsBinsAll,1)*2],mean(Stimulus2EpoksTimesDurationsBinsAll(:,3,:),3),std(Stimulus2EpoksTimesDurationsBinsAll(:,3,:),0,3)/sqrt(size(Stimulus2EpoksTimesDurationsBinsAll(:,3,:),3)),'r','LineStyle','none');
    xlabel('Time')
    ylabel('Investigation time')
    title('Long bouts along time')
    hold off; 
    
    subplot(3,3,5); 
    hold on;
    for i=1:length(TransitionTimesForHistAllForRaster)
      CurrentTransitionTimes=[];
      CurrentTransitionTimes=TransitionTimesForHistAllForRaster{i};
      RasterRowNum=[];
      RasterRowNum(1:1:length(CurrentTransitionTimes))=i;
      scatter(CurrentTransitionTimes,RasterRowNum,'b'); 
    end
    ax1 = gca;
    ax1.YColor = 'b';
    xlabel('Time');
    xlim([0 size(TransitionTimesForHistAll,2)*20*30+10*30])
    ylim([0 i+2])
    ylabel(ax1,'Animal number');
    title('Transitions alomg time');
    ax1 = gca;
    ax1_pos = ax1.Position;
    ax2 = axes('Position',ax1_pos,'XAxisLocation','bottom','YAxisLocation','right','Color','none');
    ax2.YColor = 'r';
    line(10*30:20*30:size(TransitionTimesForHistAll,2)*20*30,mean(TransitionTimesForHistAll),'Parent',ax2,'Color','r','LineWidth',2)
    xlim([0 size(TransitionTimesForHistAll,2)*20*30+10*30])
    ylabel(ax2,'average transitions number');
    hold off; 
    
    subplot(3,3,6); 
    hold on;
    bar([1:1:size(TransitionsAlongTimeInMinutes,2)],mean(TransitionsAlongTimeInMinutes,1),0.4,'b');
    errorbar([1:1:size(TransitionsAlongTimeInMinutes,2)],mean(TransitionsAlongTimeInMinutes,1),std(TransitionsAlongTimeInMinutes,0,1)/sqrt(size(TransitionsAlongTimeInMinutes,1)),'b','LineStyle','none');
    xlabel('Time');
    ylabel('Transitions number');
    title('Number of transitions along time');
    hold off; 
    
    subplot(3,3,7); 
    hold on;
    imagesc(1:size(Stimulus1InteractionHeattMap,2),1:size(Stimulus1InteractionHeattMap,1),Stimulus1InteractionHeattMap); axis tight;
    set(gca,'ydir','normal');
    xlabel('Time (frames)');
    ylabel('Animal number');
    colorscale = colorbar;
    colorscale.Label.String = 'Interaction duration in Sec (30Hz recording)';
    colorbar('Ticks',[0,3,6,9,12,15,18,20],'TickLabels',{'0','3','6','9','12','15','18','>=20'})
    caxis([0,20]);
    title(['Durations of interactions  with stimulus 1 along time (30Hz recording']);
    hold off; 

    subplot(3,3,8); 
    hold on;
    imagesc(1:size(Stimulus2InteractionHeattMap,2),1:size(Stimulus2InteractionHeattMap,1),Stimulus2InteractionHeattMap); axis tight;
    set(gca,'ydir','normal');
    xlabel('Time (frames)');
    ylabel('Animal number');
    colorscale = colorbar;
    colorscale.Label.String = 'Interaction duration in Sec (30Hz recording)';
    colorbar('Ticks',[0,3,6,9,12,15,18,20],'TickLabels',{'0','3','6','9','12','15','18','>=20'})
    caxis([0,20]);
    title(['Durations of interactions  with stimulus 2 along time (30Hz recording']);
    hold off; 
    
    TotalTimeStim1_Short_Interactions_AlongTime=squeeze(Stimulus1EpoksTimesDurationsBinsAll(1:4,1,:))';
    TotalTimeStim2_Short_Interactions_AlongTime=squeeze(Stimulus2EpoksTimesDurationsBinsAll(1:4,1,:))';
    TotalTimeStim1_Long_Interactions_AlongTime=squeeze(Stimulus1EpoksTimesDurationsBinsAll(1:4,3,:))';
    TotalTimeStim2_Long_Interactions_AlongTime=squeeze(Stimulus2EpoksTimesDurationsBinsAll(1:4,3,:))';
    TransitionsAlongTimeInMinutes_Minute1=TransitionsAlongTimeInMinutes(:,1:4);
   
end



