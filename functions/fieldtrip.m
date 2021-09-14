function [td] = fieldtrip(data,params)

td = [];
td.fsample = params.fsample;

if ~iscell(data)
    time = 1/td.fsample:1/td.fsample:(size(data,1)-1)/td.fsample;
    td.time{1} = time;     
    td.trial{1} = data';
    td.sampleinfo(1,:) = [1 length(td.time{1})];
    for cnt = 1 : size(data,2)
        td.label{cnt} = sprintf('CH%d',cnt);
    end   
    else
    if size(data{1},1)>size(data{1},2)
        time = 1/td.fsample:1/td.fsample:size(data{1},1)*1/td.fsample;
    else
        time = 1/td.fsample:1/td.fsample:size(data{1},2)*1/td.fsample; 
    end
    for cnt = 1 : length(data)
        td.time{cnt} = time;
        if size(data{cnt},2)==length(time)
            td.trial{cnt} = data{cnt};
        else
            td.trial{cnt} = data{cnt}';            
        end
        if cnt>1
            td.sampleinfo(cnt,:) = [td.sampleinfo(cnt-1,2)+1 length(td.time{1})*cnt];
        else
            td.sampleinfo(cnt,:) = [1 length(td.time{1})];            
        end
    end
    for cnt = 1 : size(td.trial{1},1)
        td.label{cnt,1} = sprintf('CH%d',cnt);
    end       
end
