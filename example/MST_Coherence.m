function [max_time, max_val, output_data, f, slice] = MST_Coherence(data1, data2, fs, overlap, scope)
    output_data = [];
    start_i = 1;
    end_i = length(data2);
    window_len = end_i;
    complete = 0;
    while complete==0
        processed1 = data1(start_i:end_i);
        [processed, w] = mscohere(processed1,data2);
        output_data = [output_data,processed];
        start_i = start_i+window_len-overlap;
        end_i = end_i+window_len-overlap;
        if end_i > length(data1)
            complete = 1;
        end
    end
    f = linspace(0, fs/2, length(w));
    f_ind1 = find(abs(f-scope(1))==min(abs(f-scope(1))));
    f_ind2 = find(abs(f-scope(2))==min(abs(f-scope(2))));
    output_data = output_data(f_ind1:f_ind2,:);
    time_axis = linspace(1/fs,(length(data1)-length(data2))/fs,length(output_data(1,:)));
    ind = find(sum(output_data,1) == max(sum(output_data,1)));
    max_time = time_axis(ind);
    max_val = max(sum(output_data,1));
    slice = output_data(:,ind);
end
