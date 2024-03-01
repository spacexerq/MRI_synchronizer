load('C:\work\test_sequence\gate_ADC_array.mat')
load('C:\work\test_sequence\gate_RF_amp.mat')
load('C:\work\test_sequence\gate_tr_sw_array.mat')

gate_array = zeros(512,5);
row_counter = 1;
cycle_counter = 1000/20;  %1us / 20ns - minimum num of cycles
% sum = 0;
for i =1:(height(gate_RF_amp_array)-1)
    if (((gate_RF_amp_array(i,2)) == (gate_RF_amp_array(i+1,2))) && ((gate_ADC_array(i,2)) == (gate_ADC_array(i+1,2))) && ((gate_tr_sw_array(i,2)) == ((gate_tr_sw_array(i+1,2)))))
        if i < height(gate_RF_amp_array)-1
            cycle_counter = cycle_counter + 1000/20;
        else
            gate_array(row_counter,1) = (gate_RF_amp_array(i,2));
            gate_array(row_counter,2) = (gate_tr_sw_array(i,2));
            gate_array(row_counter,3) = (gate_ADC_array(i,2));
            gate_array(row_counter,4) = 0;
            gate_array(row_counter,5) = (cycle_counter+(1000/20));
%             sum = sum + cycle_counter+1000/20;
            cycle_counter = 1000/20;  %1us / 20ns - minimum num of cycles
            row_counter = row_counter +1;

        end
        
    else
        gate_array(row_counter,1) = (gate_RF_amp_array(i,2));
        gate_array(row_counter,2) = (gate_tr_sw_array(i,2));
        gate_array(row_counter,3) = (gate_ADC_array(i,2));
        gate_array(row_counter,4) = 0;
        gate_array(row_counter,5) = (cycle_counter);
        row_counter = row_counter +1;
%         sum = sum + cycle_counter;
        cycle_counter = 1000/20;  %1us / 20ns - minimum num of cycles
    
        
    end %if
end %for
gate_array(row_counter:height(gate_array),:)=[]
title ={'RF_GATE' ,'SWITCH_GATE','ADC_GATE','CYCLES'};
csvwrite('gate_RF_SWITCH_ADC_GRU_CYCLES.csv',gate_array); 

% disp(sum*20);