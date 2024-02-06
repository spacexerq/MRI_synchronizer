%% Read a sequence file
% A sequence can be loaded from the open MR file format using the |read|
% method.
% seq_name='TSE_LF_16x16_TR_200_1sl.seq'; 
%seq_name='epi_rs.seq';
%seq_name='spiral.seq';
%seq_name='trufi.seq'; 
%seq_name='epi_se.seq';
%seq_name='../tests/gre.seq';
%seq_name='../tests/epi_rs.seq';

% sys = mr.opts('B0', 1.5); % we need system here if we want 'detectRFuse' to detect fat-sat pulses
% seq=mr.Sequence(sys);
% seq.read(seq_name,'detectRFuse');
%rep = seq.testReport;


% time=[];
% B1_array=[];
% for i=1:length(seq.blockDurations)
%     b1= seq.getBlock(i)
%     rf= b1.rf;
%     time1=rf.t;
%     time=[time; time1];
% %     RF=abs(rf.signal);
% %     B1_array=[B1_array; RF]
% end
%rf=b1.rf
%% Extract RF array
seq_name='SE_rfdeath_5000.seq'; 
sys = mr.opts('B0', 1.5); 
seq=mr.Sequence(sys);
seq.read(seq_name,'detectRFuse');
%In Sequence.m put the pause execution at line 1230
seq.plot()
time=round((time_for_rf+0.0000005)*1000000);%учесть сдвиг на 0.5мкс
% figure
% plot (time_for_rf, signal_array);
% plot (time_for_rf, phase_array);

RF_extracted=[];
RF_extracted= [signal_array,phase_array,time];
TR=0.200;
RF_raster_time = 0.000001;
N_rf=round(TR/RF_raster_time);

RF_extracted_new=zeros(N_rf,3);

    for k=1:length(RF_extracted)
        if RF_extracted(k,1)>0
           n=RF_extracted(k,3);
           RF_extracted_new(n,:)= RF_extracted(k,:);
           if n==N_rf
               break
           end

        end
    end

 for l=1:length(RF_extracted_new) 
        if RF_extracted_new(l,3)==0;
           RF_extracted_new(l,3)= l;% переводим в термины шагов длительностью 1 мкс
        end
 end

%Format of RF events:
%id amplitude mag_id phase_id time_shape_id delay freq phase
%..        Hz   ....     ....          ....    us   Hz   rad

RF1=[227.736 1 2 0 5000 0 0];% from .seq 
RF2=[384.179 3 4 0 5000 0 1.5708];
 %Необходимо переписать код, сделать возможным выгрузку задержек и сдвигов
 %частот сразу из .seq файла, внести столбец со сдвигом по частоте
 %заполнения



%% Plot
 figure, 
 subplot(6,1,1);
 plot (RF_extracted_new(:,3)*0.000001,RF_extracted_new(:,1))


%% Extract Grad array
seq_name='SE_rfdeath_5000.seq'; 
sys = mr.opts('B0', 1.5); % we need system here if we want 'detectRFuse' to detect fat-sat pulses
seq=mr.Sequence(sys);
seq.read(seq_name,'detectRFuse');
%In testReport.m put the pause execution at line 200
rep = seq.testReport;
%figure; plot(common_time, gw_ct');
grad_extracted=[];
grad_extracted(:,1:3)=gw_ct';
grad_extracted(:,4)=round(common_time,5);
TR=0.200;
grad_extracted(length(grad_extracted)+1,:)=[0,0,0,TR];
%% Plot
 subplot(6,1,2);
 plot(grad_extracted(:,4), grad_extracted(:,1:3)/10000);


%% gate for rf_amplif
B=load("Sequence_from_seq_SE_rfdeath_5000.mat")
gate_RF_amp=[];
block_raster_time = 0.00001;
RF_raster_time = 0.000001;
gate_RF_amp(:,1)= B.A(:,2);%*block_raster_time;
total_duration_in_blocks = sum (B.A(:,2));

% gates in term of blocks
for n=1:length(B.A(:,3))
    if B.A(n,3)>0
    gate_RF_amp(n,2) = 1;
    else
    gate_RF_amp(n,2) = 0;
    end
end    

scale=round(block_raster_time/RF_raster_time);

% gates in term of RF raster
gate_RF_amp_array=[];
for n=1:length(B.A(:,3))
    if n==1
        for i=1:gate_RF_amp(n,1)*scale;
            gate_RF_amp_array(i,1)= i;
            gate_RF_amp_array(i,2)= gate_RF_amp(n,2);
        end
        count=gate_RF_amp(n,1)*scale;
    else
        for i=count+1:count+gate_RF_amp(n,1)*scale
            gate_RF_amp_array(i,1)= i;
            gate_RF_amp_array(i,2)= gate_RF_amp(n,2);
        end
        count=count+gate_RF_amp(n,1)*scale;
    end
end
%% Plot
subplot(6,1,3);

plot(gate_RF_amp_array(:,1)/1000000,gate_RF_amp_array(:,2)*400)

%% gate for ADC
B=load("Sequence_from_seq_SE_rfdeath_5000.mat")
gate_ADC=[];
block_raster_time = 0.00001;
gate_ADC(:,1)= B.A(:,2);%*block_raster_time;
%total_duration_blocks = sum (B.A(:,2));
TR=0.200;

%gates in terms of blocks
for n=1:length(B.A(:,7))
    if B.A(n,7)>0
    gate_ADC(n,2) = 1;
    else
    gate_ADC(n,2) = 0;
    end
end 

% Correction for grad rize and fall times - надо включить АЦП только во
% время полки градиента!
%Руками из seq файла взяла freq_enc_grad
freq_enc_grad=[ 8  45454.5  10 2200  10   0];
freq_enc_grad_rize_time_in_blocks=freq_enc_grad(1,3)/1000000/block_raster_time;%кол-во блоков на рост
freq_enc_grad_fall_time_in_blocks=freq_enc_grad(1,5)/1000000/block_raster_time;%кол-во блоков на спад

% n_ADC=sum(gate_ADC(:,2) == 1); %number of ADC events
% new_number_of_blocks=length(gate_ADC(:,1))+n_ADC*2;

gate_ADC_corr=[];
for n=1:length(gate_ADC)
    if gate_ADC(n,2)==0
    gate_ADC_corr= [gate_ADC_corr;gate_ADC(n,:)]
    else
    gate_ADC(n,1)=gate_ADC(n,1)-freq_enc_grad_rize_time_in_blocks-freq_enc_grad_fall_time_in_blocks;
    gate_ADC_corr= [gate_ADC_corr;[freq_enc_grad_rize_time_in_blocks,0];gate_ADC(n,:);[freq_enc_grad_fall_time_in_blocks,0]]
    end    
end
gate_ADC= gate_ADC_corr;

% gates in term of RF raster
scale=round(block_raster_time/RF_raster_time);%

gate_ADC_array=[];
for n=1:length(gate_ADC(:,1))
    if n==1
        for i=1:gate_ADC(n,1)*scale;
            gate_ADC_array(i,1)= i;
            gate_ADC_array(i,2)= gate_ADC(n,2);
        end
        count=gate_ADC(n,1)*scale;
    else
        for i=count+1:count+gate_ADC(n,1)*scale
            gate_ADC_array(i,1)= i;
            gate_ADC_array(i,2)= gate_ADC(n,2);
        end
        count=count+gate_ADC(n,1)*scale;
    end
end


%% Plot

subplot(6,1,4);
plot(gate_ADC_array(:,1)/1000000,gate_ADC_array(:,2)*400)

%% gate for tr switch
%Переделываем из гейта на АЦП
B=load("Sequence_from_seq_SE_rfdeath_5000.mat")
gate_ADC=[];
block_raster_time = 0.00001;
gate_ADC(:,1)= B.A(:,2);%*block_raster_time;
%total_duration_blocks = sum (B.A(:,2));
TR=0.200;

%gates in terms of blocks
for n=1:length(B.A(:,7))
    if B.A(n,7)>0
    gate_ADC(n,2) = 1;
    else
    gate_ADC(n,2) = 0;
    end
end 

% Correction for grad rize and fall times
%Руками из seq файла взяла freq_enc_grad
freq_enc_grad=[ 8  45454.5  10 2200  10   0];
freq_enc_grad_rize_time_in_blocks=freq_enc_grad(1,3)/1000000/block_raster_time;%кол-во блоков на рост
freq_enc_grad_fall_time_in_blocks=freq_enc_grad(1,5)/1000000/block_raster_time;%кол-во блоков на спад

% n_ADC=sum(gate_ADC(:,2) == 1); %number of ADC events
% new_number_of_blocks=length(gate_ADC(:,1))+n_ADC*2;

gate_ADC_corr=[];
for n=1:length(gate_ADC)
    if gate_ADC(n,2)==0
    gate_ADC_corr= [gate_ADC_corr;gate_ADC(n,:)]
    else
    gate_ADC(n,1)=gate_ADC(n,1)-freq_enc_grad_rize_time_in_blocks-freq_enc_grad_fall_time_in_blocks;
    gate_ADC_corr= [gate_ADC_corr;[freq_enc_grad_rize_time_in_blocks,0];gate_ADC(n,:);[freq_enc_grad_fall_time_in_blocks,0]]
    end    
end
gate_ADC= gate_ADC_corr;

gate_tr_sw=gate_ADC;
scale=round(block_raster_time/RF_raster_time);%можно и ADC raster поставить, но надо ли?

%Учет задержки до АЦП, тк свич может звенеть 100 мкс????
tr_sw_fall=0.0001;%s 100 микросекунд - добавть правильную!!!!!
tr_sw_fall_in_blocks=tr_sw_fall/block_raster_time;%кол-во блоков на задержку

gate_tr_sw_corr=[];
temp=zeros(1,2); 
for n=1:length(gate_tr_sw)
    if gate_tr_sw(n,2)==0
    temp(1,1)=temp(1,1)+gate_tr_sw(n,1); %сложила все нулевые блоки, надо так сделать и для остальных гейтов
    %gate_tr_sw_corr= [gate_tr_sw_corr;temp]
    else
    temp(1,1)=temp(1,1)-tr_sw_fall_in_blocks;
    gate_tr_sw_with_fall=gate_tr_sw(n,:);
    gate_tr_sw_with_fall(1,1)=gate_tr_sw_with_fall(1,1)+tr_sw_fall_in_blocks;
    gate_tr_sw_corr= [gate_tr_sw_corr;temp;gate_tr_sw_with_fall];
    temp=zeros(1,2); 
    end    
    if n==length(gate_tr_sw) && temp(1,1)>0
    gate_tr_sw_corr= [gate_tr_sw_corr;temp];% Чтобы последняя серия нулевых блоков тоже шла в последовательность гейтов
    end
end
gate_tr_sw=gate_tr_sw_corr;

% 
% % gates in term of RF raster
gate_tr_sw_array=[];
for n=1:length(gate_tr_sw(:,1))
    if n==1
        for i=1:gate_tr_sw(n,1)*scale;
            gate_tr_sw_array(i,1)= i;
            gate_tr_sw_array(i,2)= gate_tr_sw(n,2);
        end
        count=gate_tr_sw(n,1)*scale;
    else
        for i=count+1:count+gate_tr_sw(n,1)*scale
            gate_tr_sw_array(i,1)= i;
            gate_tr_sw_array(i,2)= gate_tr_sw(n,2);
        end
        count=count+gate_tr_sw(n,1)*scale;
    end
end

%% Plot
%%
figure,
subplot(3,1,1);

