
function [data, digitalChannels, time] = readIgor_withDigital2()
%f=fopen(filename);
%cd('C:\data');
[name, path] = uigetfile('*.*');
cd(path)
f = fopen(name);
fprintf('file opened: %s \n', name);

%f=fopen('C:\Documents and Settings\hestrin\Desktop\081210_A.004')
stringsize=fread(f,1,'single');
wavesize=fread(f,1,'single');
headers = transpose(char(fread(f,stringsize,'char*1')));
headers = strrep(headers, ';', ' ');
hwave=fread(f,wavesize,'single');
data=fread(f,inf,'integer*2', 'l');
dataSize = size(data);

%extract necessary parameters of data acquisition 
trace_num = sscanf(headers, ['acquired:' '%d']);

samplesPerTrace_index = strfind(headers, 'samples:');
samplesPerTrace = sscanf(headers(samplesPerTrace_index:end), ['samples:' '%d']);

frequency_index = strfind(headers, 'freq:');
frequency = sscanf(headers(frequency_index:end), ['freq:' '%d']); 

total_chan_num_index = strfind(headers, 'total_chan_num:');                                 % index where the total_chan_num first occurs in the header string
total_chan_num = sscanf(headers(total_chan_num_index:end), ['total_chan_num:' '%d']);       % store the actual channel number in this variable: THIS INCLUDES DUMMY CHANNELS

adc_status0_index = strfind(headers, 'adc_status0>');
adc_flag(1) = sscanf(headers(adc_status0_index:end), ['adc_status0>' '%d']); 

adc_status1_index = strfind(headers, 'adc_status1>');
adc_flag(2) = sscanf(headers(adc_status1_index:end), ['adc_status1>' '%d']); 

adc_status2_index = strfind(headers, 'adc_status2>');
adc_flag(3) = sscanf(headers(adc_status2_index:end), ['adc_status2>' '%d']); 

adc_status3_index = strfind(headers, 'adc_status3>');
adc_flag(4) = sscanf(headers(adc_status3_index:end), ['adc_status3>' '%d']); 

adc_status4_index = strfind(headers, 'adc_status4>');
adc_flag(5) = sscanf(headers(adc_status4_index:end), ['adc_status4>' '%d']); 

adc_status5_index = strfind(headers, 'adc_status5>');
adc_flag(6) = sscanf(headers(adc_status5_index:end), ['adc_status5>' '%d']); 

adc_status6_index = strfind(headers, 'adc_status6>');
adc_flag(7) = sscanf(headers(adc_status6_index:end), ['adc_status6>' '%d']); 

adc_status7_index = strfind(headers, 'adc_status7>');
adc_flag(8) = sscanf(headers(adc_status7_index:end), ['adc_status7>' '%d']); 

adcs = ['adc_replace0:'; 'adc_replace1:'; 'adc_replace2:'; 'adc_replace3:'; 'adc_replace4:'; 'adc_replace5:'; 'adc_replace6:'; 'adc_replace7:'];
active_adcs = adcs(adc_flag==1, :);
channel_num = sum(adc_flag);                                                                % THIS INCLUDES ONLY NON-DUMMY CHANNELS

gain = zeros(channel_num, 1);
for i=1:channel_num
    adc = active_adcs(i, :);
    adc_gain_str = strrep(adc, 'replace', 'gain');
    gain_index = strfind(headers, adc_gain_str);
    gain(i, 1) = 3.2*sscanf(headers(gain_index:end), [adc_gain_str '%d']);
end

% if data wasn't saved properly and acquired reads as 0, figure out how
% many traces are there
if trace_num == 0 
    trace_num = (dataSize(1)/total_chan_num)/samplesPerTrace;
end


data = reshape(data, total_chan_num, dataSize(1)/total_chan_num);

data = data(1:channel_num+1, :);
digiData = squeeze(data(channel_num+1, :));
digiData(digiData < 0) = digiData(digiData<0) + 2*2^15;
data(channel_num + 1, :) = digiData;
gain = [gain; 1];

gain = repmat(gain, 1, dataSize(1)/total_chan_num);
data = data./gain;
    
[channels samples] = size(data);
traceLength = samples/trace_num;
data = reshape(data, channels, traceLength, trace_num);

time = 0: 1/frequency: (traceLength-1)/frequency;


digidata = squeeze(data(end, :, :));
ttls = dec2bin(digidata);
ttls=double(ttls-'0');
digitalChannels = zeros(size(ttls));
for i = 1:size(ttls, 2)
    digitalChannels(:, i) = ttls(:, size(ttls, 2) - i + 1);
end

digitalChannels = reshape(digitalChannels', size(ttls, 2), traceLength, trace_num);

