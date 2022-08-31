classdef Expt3_Keypad_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        Result                          matlab.ui.control.TextArea
        TabGroup                        matlab.ui.container.TabGroup
        Parameters                      matlab.ui.container.Tab
        YoucanchangetheParametersbeforedecodingPanel  matlab.ui.container.Panel
        GridLayout2                     matlab.ui.container.GridLayout
        FilterLengthEditField           matlab.ui.control.NumericEditField
        FilterLengthEditFieldLabel      matlab.ui.control.Label
        NoiseAmplitudeEditField         matlab.ui.control.NumericEditField
        NoiseAmplitudeEditFieldLabel    matlab.ui.control.Label
        TextArea_2                      matlab.ui.control.TextArea
        NEXTButton                      matlab.ui.control.Button
        SamplePointsEditField           matlab.ui.control.NumericEditField
        SamplePointsEditFieldLabel      matlab.ui.control.Label
        TimeDurationofeachbitEditField  matlab.ui.control.NumericEditField
        TimeDurationofeachbitEditFieldLabel  matlab.ui.control.Label
        SamplingFrequencyinkHzEditField  matlab.ui.control.NumericEditField
        SamplingFrequencyinkHzEditFieldLabel  matlab.ui.control.Label
        KeyPad_Tab                      matlab.ui.container.Tab
        GridLayout                      matlab.ui.container.GridLayout
        Button_1                        matlab.ui.control.Button
        BackButton                      matlab.ui.control.Button
        DButton                         matlab.ui.control.Button
        CButton                         matlab.ui.control.Button
        BButton                         matlab.ui.control.Button
        AButton                         matlab.ui.control.Button
        EnterButton                     matlab.ui.control.Button
        ClearButton                     matlab.ui.control.Button
        KEYPADLabel                     matlab.ui.control.Label
        Button_12                       matlab.ui.control.Button
        Button_11                       matlab.ui.control.Button
        Button_10                       matlab.ui.control.Button
        Button_9                        matlab.ui.control.Button
        Button_8                        matlab.ui.control.Button
        Button_7                        matlab.ui.control.Button
        Button_6                        matlab.ui.control.Button
        Button_5                        matlab.ui.control.Button
        Button_4                        matlab.ui.control.Button
        TextArea                        matlab.ui.control.TextArea
        Button_3                        matlab.ui.control.Button
        Button_2                        matlab.ui.control.Button
        EncodedGraphsTab                matlab.ui.container.Tab
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
        ContextMenu                     matlab.ui.container.ContextMenu
        ClearstheentirescreenanderasestheprocessMenu  matlab.ui.container.Menu
        ContextMenu2                    matlab.ui.container.ContextMenu
        BackspacedeletesthedatacorrespondingtolastinputkeyMenu  matlab.ui.container.Menu
        ContextMenu3                    matlab.ui.container.ContextMenu
        Menu                            matlab.ui.container.Menu
        ContextMenu4                    matlab.ui.container.ContextMenu
        FixtheabovevaluesandstartdecodingMenu  matlab.ui.container.Menu
    end

    
    properties (Access = public)
        String = "";% This string stores the values entered sequentially
%         rowvec = [];
%         colvec = [];
        Fs = 5e3; %sampling frequency
        L = 0.1 ;   %Time Duration of individual bit in seconds
        N = 512 ;  %Sample Points
        i = 0;
        x_t = [];
        filt_len = 40;
        noise_amp = 0.5;
%         ftime = figure("Name", "Encoded signal in time domain");
%         ffreq = figure("Name", "Encoded signal in frequency domain");
    end
    
    methods (Access = public)
        
        function results = display_res(app)
            disp("Entered input: "+app.String);
%             disp("Row Vector: ");
%             disp(app.rowvec);
            
            results = app.String;
            encoded = app.String;
            decoded = app.Decode();
            if(encoded==decoded)
                disp("Successfully Decoded:");
                disp(char(9)+"KeyWord Encoded = "+encoded);
                disp(char(9)+"KeyWord Decoded = "+decoded);
                app.Result.Value = "Decoding Successful: "+newline+"Keyword decoded = "+decoded;
                app.Result.Visible = "on";
            else
                disp("Decoding Unsuccessful:");
                disp(char(9)+"KeyWord Encoded = "+encoded);
                disp(char(9)+"KeyWord Decoded = "+decoded);
                app.Result.Value = "Decoding Unsuccessful: "+newline+"Keyword decoded = "+decoded;
                app.Result.Visible = "on";
            end
        end
        
        function h_n = BPFilt(app, f_c, Fs)
            %designing bandpass filter bank
            w_c = f_c*(pi/(Fs/2));  %Digitizing central frequency
            wid = app.filt_len; %   Length of filter
            h_n = zeros(1, wid+1);
            for n=1:(wid)
               h_n(n) = cos(w_c*(n-1)); 
            end
            omega = linspace(-pi, pi, 750); %sweep of omega
            H_f = freqz(h_n, 1, omega); % Frequency response of the filter
            beta = 1/max(abs(H_f)); %Choosing beta to make the maximum value of frequency response magnitude as 1
            h_n = h_n*beta;
        end
    end
    
    methods (Access = private)
        
        function Encode(app, f1, f2)
                app.Result.Visible = "off";
                
                t = ((app.i-1)*app.L):(1/app.Fs):((app.i*app.L)-(1/app.Fs));
                x1_t = sin(2*pi*f1*t) + sin(2*pi*f2*t);
                app.x_t = cat(2, app.x_t, x1_t);  %Input signal
                
                N1 = app.N*app.i;
                X_w = fftshift(fft(app.x_t, N1))/N1; %fft sprectrum shifting and normalization
                f = linspace(-1/2,1/2,N1)*(app.Fs);
                
                ax1 = app.UIAxes;
                plot(ax1, (0:uint64(app.i*app.L*app.Fs)-1), app.x_t(1:uint64(app.i*app.L*app.Fs)), "LineWidth", 1.25);               
                title(ax1,"Encoded Signal in Time Domain", "fontsize", 16, "Color", 'r');
                axis(ax1, "tight");
                grid(ax1, "on");
                grid(ax1, "minor");
                xlabel(ax1,'time (s)', "fontsize", 12);
                ylabel(ax1,'Amplitude', "fontsize", 12);
                legend(ax1,"Input Signal");
                
                ax2 = app.UIAxes2;
                plot(ax2, f, abs(X_w), "LineWidth", 1.25);
                xlabel(ax2, 'Frequency (in Hz)', "fontsize", 12);
                ylabel(ax2, '|X(f)/N|', "fontsize", 12);
                title(ax2, "Encoded Signal Frequency spectrum","fontsize", 16, "Color", 'r');
                legend(ax2, "Input Signal");
                axis(ax1, "tight");
                grid(ax1, "on");
                grid(ax1, "minor");
                
        end
        
        
        function rslt = Decode(app)
                freqvec = [697, 770, 852, 941, 1209, 1336, 1477, 1633];
                rslt = "";
                len = app.L;
                num_bit = length(app.x_t)/(len*app.Fs); %Number of bits in the input
                t = 0:1/app.Fs:(num_bit*len-1/app.Fs);    %len is individual bit length
                
                N1 = app.N*num_bit;
                
                f = linspace(-1/2,1/2,N1)*(app.Fs);
                
                noise = app.noise_amp*randn(1,length(app.x_t)); %%white noise
                x_noise = app.x_t + noise; %%noisy input
                
                %Time domain of Input signal
                figure("Name", "Plotting Received signals");
                subplot(2,2,1);
                plot(t, app.x_t, "LineWidth", 1.25);
                xlabel('time (in seconds)', "fontsize", 12);
                ylabel('x(t)', "fontsize", 12);
                title("Time spectrum of received signal", "fontsize", 14, "Color", 'r');
                axis tight;
                grid on;
                grid minor;
                
                %Time domain of Noise Input signal
                subplot(2,2,3);
                plot(t, x_noise, "LineWidth", 1.25);
                xlabel('time (in seconds)', "fontsize", 12);
                ylabel('x(t)', "fontsize", 12);
                title("Time spectrum of received noisy signal", "fontsize", 14, "Color", 'r');
                axis tight;
                grid on;
                grid minor;
                
                
                %Freq Spectrum of Input signal
                subplot(2,2,2);
                plot(f, abs(fftshift(fft(app.x_t, N1))/N1), "LineWidth", 1.25);
                xlabel('Frequency (in Hz)', "fontsize", 12);
                ylabel('|X(f)/N|', "fontsize", 12);
                title("Frequency spectrum of received signal", "fontsize", 14, "Color", 'r');
                axis tight;
                grid on;
                grid minor;
                
                %Freq Spectrum of noisy Input signal
                subplot(2,2,4);
                plot(f, abs(fftshift(fft(x_noise, N1))/N1), "LineWidth", 1.25);
                xlabel('Frequency (in Hz)', "fontsize", 10);
                ylabel('|X(f)/N|', "fontsize", 10);
                title("Frequency spectrum of received noise added signal", "fontsize", 14, "Color", 'r');
                axis tight;
                grid on;
                grid minor;
                
                
                for k = 1:num_bit
                    maxmagval = zeros(1, length(freqvec));
                    maxmagval_n = zeros(1, length(freqvec));
                    res = "";
                    x1_t = app.x_t(uint64((k-1)*len*app.Fs+1):uint64(k*len*app.Fs));
                    x1_noise = x_noise(uint64((k-1)*len*app.Fs+1):uint64(k*len*app.Fs));
                    t1 = ((k-1)*len):1/app.Fs:(k*len-1/app.Fs);
                    f_1 = linspace(-1/2,1/2,app.N)*app.Fs;
                    
                    %Time domain of Input signal
                    f_time = figure("Name", "Bit "+int2str(k)+" Time Domain");
                    subplot(3,3,1);
                    plot(t1, x1_t, "LineWidth", 1.25);
                    xlabel('time (in seconds)', "fontsize", 9);
                    ylabel('x(t)', "fontsize", 9);
                    title("Bit "+int2str(k)+" Time Spectrum", "fontsize", 12.5, "Color", 'r');
                    legend("input signal");
                    axis tight;
                    grid on;
                    grid minor;
            
                    %Time domain of Noise Input signal
                    f_time_noise = figure("Name", "Bit "+int2str(k)+" Time Domain Noisy");
                    subplot(3,3,1);
                    plot(t1, x1_noise, "LineWidth", 1.25);
                    xlabel('time (in seconds)', "fontsize", 9);
                    ylabel('x(t)', "fontsize", 9);
                    title("Bit "+int2str(k)+" Time Spectrum Noise Added", "fontsize", 11.5, "Color", 'r');
                    legend("input signal");
                    axis tight;
                    grid on;
                    grid minor;
            
            
                    %Freq Spectrum of Input signal
                    f_freq = figure("Name", "Bit "+int2str(k)+" Freq Domain");
                    subplot(3,3,1);
                    plot(f_1, abs(fftshift(fft(x1_t, app.N))/app.N), "LineWidth", 1.25);
                    xlabel('Frequency (in Hz)', "fontsize", 9);
                    ylabel('|X(f)/N|', "fontsize", 9);
                    title("Bit "+int2str(k)+" Frequency Spectrum", "fontsize", 11, "Color", 'r');
                    legend("input signal");
                    axis tight;
                    grid on;
                    grid minor;
            
                    %Freq Spectrum of noisy Input signal
                    f_freq_noise = figure("Name", "Bit "+int2str(k)+" Freq Domain Noisy");
                    subplot(3,3,1);
                    plot(f_1, abs(fftshift(fft(x1_noise, app.N))/app.N), "LineWidth", 1.25);
                    xlabel('Frequency (in Hz)', "fontsize", 9);
                    ylabel('|X(f)/N|', "fontsize", 9);
                    title("Bit "+int2str(k)+" Frequency Spectrum Noise Added", "fontsize", 10.5, "Color", 'r');
                    legend("input signal");
                    axis tight;
                    grid on;
                    grid minor;
                    
                    %         f_freq_H = figure("Name", "Filter Response");
                    
                    for j=1:length(freqvec)
                %         fpass = [freqvec(i)-10, freqvec(i)+10];   %for using inbuilt filter
                %         y_n = bandpass(x_t, fpass, Fs);
                        h_n = app.BPFilt(freqvec(j), app.Fs);
            
                        y_n = filtfilt(h_n, 1, x1_t);    %Filtering the signal using filtfilt command
            
                        y_noise = filtfilt(h_n, 1, x1_noise);    %Filtering the noisy signal using filtfilt command
            
                        Y_w = fftshift(fft(y_n, app.N))/app.N; %fft sprectrum shifting and normalization
            
                        Y_w_noise = fftshift(fft(y_noise, app.N))/app.N; %fft sprectrum shifting and normalization
            
            %             %         Plotting the frequency domain representation of the h_n
            %             figure(f_freq_H);
            %             subplot(3,3,j);
            %             omega = linspace(-pi, pi, 750); %sweep of omega
            %             H_f = freqz(h_n, 1, omega); % Frequency response of the filter
            %             
            %             plot(omega/pi*(Fs/2), 20*log10(abs(H_f)), "LineWidth", 1.3);
            %             xlabel("Frequency (Hz)", "fontsize", 10);
            %             ylabel('Magnitude (dB)',"fontsize", 10);
            %             title("Frequency Response of BPF with fc = "+freqvec(j)+" Hz", "fontsize", 14, "Color", 'r');
            %             axis tight;
            %             grid on;
            %             grid minor;
            
                        %Time Spectrum of Filtered signal
                        figure(f_time);
                        subplot(3,3, j+1);
                        plot(t1, y_n, "LineWidth", 1.25);   %Plotting for first 128 points
                        xlabel('Time (in seconds)', "fontsize", 9);
                        ylabel('y(t)', "fontsize", 9);
                        title("Time spectrum of BPF filtered", "fontsize", 11.5, "Color", 'r');
                        legend("f_c = "+int2str(freqvec(j))+ " Hz");
                        axis tight;
                        grid on;
                        grid minor;
            
                        %Time Spectrum of Noisy Filtered signal
                        figure(f_time_noise);
                        subplot(3,3, j+1);
                        plot(t1, y_noise, "LineWidth", 1.25);   %Plotting for 2 sampling periods
                        xlabel('Time (in seconds)', "fontsize", 9);
                        ylabel('y(t)', "fontsize", 9);
                        title("Time spectrum of BPF filtered noisy signal", "fontsize", 10, "Color", 'r');
                        legend("f_c = "+int2str(freqvec(j))+ " Hz");
                        axis tight;
                        grid on;
                        grid minor;        
            
                        %Freq Spectrum of Filtered signal
                        figure(f_freq);
                        subplot(3,3, j+1);
                        plot(f_1, abs(Y_w), "LineWidth", 1.25);
                        xlabel('Frequency (in Hz)', "fontsize", 9);
                        ylabel('|X(f)/N|', "fontsize", 9);
                        title("Frequency spectrum of BPF filtered", "fontsize", 11.5, "Color", 'r');
                        legend("f_c = "+int2str(freqvec(j))+ " Hz");
                        axis tight;
                        grid on;
                        grid minor;
            
                        %Freq Spectrum of Noisy Filtered signal
                        figure(f_freq_noise);
                        subplot(3,3, j+1);
                        plot(f_1, abs(Y_w_noise), "LineWidth", 1.25);
                        xlabel('Frequency (in Hz)', "fontsize", 9);
                        ylabel('|X(f)/N|', "fontsize", 9);
                        title("Frequency spectrum of BPF filtered noisy", "fontsize", 10, "Color", 'r');
                        legend("f_c = "+int2str(freqvec(j))+ " Hz");
                        axis tight;
                        grid on;
                        grid minor;        
            
                        [temp, ind] = max(Y_w);%finding the max freq value and storing its corresponding freq value
                        maxmagval(j) = temp;
            
                        temp_n = rms(Y_w_noise);%finding the max freq value and storing its corresponding freq value
                        maxmagval_n(j) = temp_n;
                %         disp("for h_"+int2str(i)+" filtered, max val freq is at: "+int2str(abs(f(ind))));
                    end
                
                    [tmp, I] = maxk(maxmagval, 2); %Finding the max 2 of all the corresponding values and storing the corresponding indices and from that indices finding the corresponding frequecies from freqvec
                    maxim = [freqvec(I(1)), freqvec(I(2))];
                    res1 = app.findkey(maxim);
            
                    %For noisy frequency response
                    [tmp_n, I_n] = maxk(maxmagval_n, 2); %Finding the max 2 of all the corresponding values and storing the corresponding indices and from that indices finding the corresponding frequecies from freqvec
                    maxim_n = [freqvec(I_n(1)), freqvec(I_n(2))];
                    res2 = app.findkey(maxim_n);
            
                    if(res1==res2)  
                        res = res1;
                    else 
                        res = 'W'; %Returns Incorrect result
                    end
                    rslt = rslt + res;  
                end
        end
        
        function res = findkey(app, maxim)
            if(ismember(697, maxim) && ismember(1209, maxim)) res = '1';
            elseif (ismember(697, maxim) && ismember(1336, maxim)) res = '2';
            elseif (ismember(697, maxim) && ismember(1477, maxim)) res = '3';
            elseif (ismember(697, maxim) && ismember(1633, maxim)) res = 'A'; 
            elseif (ismember(770, maxim) && ismember(1209, maxim)) res = '4';
            elseif (ismember(770, maxim) && ismember(1336, maxim)) res = '5';
            elseif (ismember(770, maxim) && ismember(1477, maxim)) res = '6';
            elseif (ismember(770, maxim) && ismember(1633, maxim)) res = 'B';%Freq Spectrum of Filtered signal
            elseif (ismember(852, maxim) && ismember(1209, maxim)) res = '7';
            elseif (ismember(852, maxim) && ismember(1336, maxim)) res = '8';
            elseif (ismember(852, maxim) && ismember(1477, maxim)) res = '9';
            elseif (ismember(852, maxim) && ismember(1633, maxim)) res = 'C';
            elseif (ismember(941, maxim) && ismember(1209, maxim)) res = '*';
            elseif (ismember(941, maxim) && ismember(1336, maxim)) res = '0';
            elseif (ismember(941, maxim) && ismember(1477, maxim)) res = '#';
            elseif (ismember(941, maxim) && ismember(1633, maxim)) res = 'D';
            end              
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: Button_1
        function Button_1Pushed(app, event)
            app.String = app.String + num2str(1); % Set
%             app.rowvec = [app.rowvec, 697];
%             app.TextArea.Value = app.String+" "+num2str(app.Fs)+" "+num2str(app.N)+" "+sprintf('%.2f', app.L); % Set
            app.TextArea.Value = app.String;
            app.i = app.i+1;
            app.Encode(697, 1209);
        end

        % Button pushed function: Button_2
        function Button_2Pushed(app, event)
            app.String = app.String + num2str(2); % Set
%             app.rowvec = [app.rowvec, 697];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(697, 1336);
        end

        % Button pushed function: Button_3
        function Button_3Pushed(app, event)
            app.String = app.String + num2str(3); % Set
%             app.rowvec = [app.rowvec, 697];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(697, 1477);
        end

        % Button pushed function: AButton
        function AButtonPushed(app, event)
            app.String = app.String + "A"; % Set
%             app.rowvec = [app.rowvec, 697];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(697, 1633);
        end

        % Button pushed function: Button_4
        function Button_4Pushed(app, event)
            app.String = app.String + num2str(4); % Set
%             app.rowvec = [app.rowvec, 770];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(770, 1209);
        end

        % Button pushed function: Button_5
        function Button_5Pushed(app, event)
            app.String = app.String + num2str(5); % Set
%             app.rowvec = [app.rowvec, 770];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(770, 1336);
        end

        % Button pushed function: Button_6
        function Button_6Pushed(app, event)
            app.String = app.String + num2str(6); % Set
%             app.rowvec = [app.rowvec, 770];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(770, 1477);
        end

        % Button pushed function: BButton
        function BButtonPushed(app, event)
            app.String = app.String + "B"; % Set
%             app.rowvec = [app.rowvec, 770];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(770, 1633);
        end

        % Button pushed function: Button_7
        function Button_7Pushed(app, event)
            app.String = app.String + num2str(7); % Set
%             app.rowvec = [app.rowvec, 852];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(852, 1209);
        end

        % Button pushed function: Button_8
        function Button_8Pushed(app, event)
            app.String = app.String + num2str(8); % Set
%             app.rowvec = [app.rowvec, 852];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(852, 1336);
        end

        % Button pushed function: Button_9
        function Button_9Pushed(app, event)
            app.String = app.String + num2str(9); % Set
%             app.rowvec = [app.rowvec, 852];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(852, 1477);
        end

        % Button pushed function: CButton
        function CButtonPushed(app, event)
            app.String = app.String + "C"; % Set
%             app.rowvec = [app.rowvec, 852];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(852, 1633);
        end

        % Button pushed function: Button_10
        function Button_10Pushed(app, event)
            app.String = app.String + "*"; % Set
%             app.rowvec = [app.rowvec, 941];
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(941, 1209);
        end

        % Button pushed function: Button_11
        function Button_11Pushed(app, event)
            app.String = app.String + num2str(0); % Set
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(941, 1336);
        end

        % Button pushed function: Button_12
        function Button_12Pushed(app, event)
            app.String = app.String + "#"; % Set
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(941, 1477);
        end

        % Button pushed function: DButton
        function DButtonPushed(app, event)
            app.String = app.String + "D"; % Set
            app.TextArea.Value = app.String; % Set
            app.i = app.i+1;
            app.Encode(941, 1633);
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            app.String = ""; %RESETTING the whole thing
            app.i = 0;
            app.TextArea.Value = app.String; % Set
            app.Fs = 5e3;
            app.L = 0.1;
            app.N = 512;
            
            app.x_t = [];
            app.filt_len = 40;
            app.noise_amp = 0.5;
            
            app.SamplingFrequencyinkHzEditField.Editable = "ON";
            app.TimeDurationofeachbitEditField.Editable = "ON";
            app.SamplePointsEditField.Editable = "ON";
            app.FilterLengthEditField.Editable = "ON";
        end

        % Button pushed function: NEXTButton
        function NEXTButtonPushed(app, event)
            app.Fs = app.SamplingFrequencyinkHzEditField.Value*1e3;
            app.L = app.TimeDurationofeachbitEditField.Value;
            app.N = app.SamplePointsEditField.Value;
            app.filt_len = app.FilterLengthEditField.Value;
            app.noise_amp = app.NoiseAmplitudeEditField.Value;
            app.x_t = [];
            app.i = 0;
            app.SamplingFrequencyinkHzEditField.Editable = "OFF";
            app.TimeDurationofeachbitEditField.Editable = "OFF";
            app.SamplePointsEditField.Editable = "OFF";
            app.NoiseAmplitudeEditField.Editable = "OFF";
            app.FilterLengthEditField.Editable = "OFF";
        end

        % Button pushed function: BackButton
        function BackButtonPushed(app, event)
            if(app.String~="")
                C = char(app.String);       %Convert string to Char
                C = C(1:length(C)-1);   %Delete the last character
                app.String = convertCharsToStrings(C);  %Convert back Char to string                
                %   Remove the previous last two given signals and x_t and
                %   encode the 2nd last signal again
                app.x_t = app.x_t(1:(app.i-1)*app.L*app.Fs);
                app.i = app.i-1;
                
                %Plotting back the previous 
                N1 = app.N*app.i;
                X_w = fftshift(fft(app.x_t, N1))/N1; %fft sprectrum shifting and normalization
                f = linspace(-1/2,1/2,N1)*(app.Fs);
                
                ax1 = app.UIAxes;
                plot(ax1, (0:(app.i*app.L*app.Fs)-1), app.x_t(1:(app.i*app.L*app.Fs)), "LineWidth", 1.25);               
                title(ax1,"Time Domain Input Signal", "fontsize", 16, "Color", 'r');
                axis(ax1, "tight");
                grid(ax1, "on");
                grid(ax1, "minor");
                xlabel(ax1,'time (s)', "fontsize", 12);
                ylabel(ax1,'Amplitude', "fontsize", 12);
                legend(ax1,"Input Signal");
                
                ax2 = app.UIAxes2;
                plot(ax2, f, abs(X_w), "LineWidth", 1.25);
                xlabel(ax2, 'Frequency (in Hz)', "fontsize", 12);
                ylabel(ax2, '|X(f)/N|', "fontsize", 12);
                title(ax2, "Encoded Frequency spectrum of input signal","fontsize", 16, "Color", 'r');
                legend(ax2, "Input Signal");
                axis(ax1, "tight");
                grid(ax1, "on");
                grid(ax1, "minor");                
            end
            app.TextArea.Value = app.String; % Set            
        end

        % Button pushed function: EnterButton
        function EnterButtonPushed(app, event)
            app.TextArea.Value = "Encoded " + app.String; % Set
            app.display_res();
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 760 662];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [143 1 561 554];

            % Create Parameters
            app.Parameters = uitab(app.TabGroup);
            app.Parameters.Title = 'Parameters';

            % Create YoucanchangetheParametersbeforedecodingPanel
            app.YoucanchangetheParametersbeforedecodingPanel = uipanel(app.Parameters);
            app.YoucanchangetheParametersbeforedecodingPanel.TitlePosition = 'centertop';
            app.YoucanchangetheParametersbeforedecodingPanel.Title = 'Signal Parameters';
            app.YoucanchangetheParametersbeforedecodingPanel.FontName = 'Arial';
            app.YoucanchangetheParametersbeforedecodingPanel.FontAngle = 'italic';
            app.YoucanchangetheParametersbeforedecodingPanel.FontWeight = 'bold';
            app.YoucanchangetheParametersbeforedecodingPanel.FontSize = 15;
            app.YoucanchangetheParametersbeforedecodingPanel.Position = [76 2 422 521];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.YoucanchangetheParametersbeforedecodingPanel);
            app.GridLayout2.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create SamplingFrequencyinkHzEditFieldLabel
            app.SamplingFrequencyinkHzEditFieldLabel = uilabel(app.GridLayout2);
            app.SamplingFrequencyinkHzEditFieldLabel.HorizontalAlignment = 'center';
            app.SamplingFrequencyinkHzEditFieldLabel.WordWrap = 'on';
            app.SamplingFrequencyinkHzEditFieldLabel.FontName = 'Calibri';
            app.SamplingFrequencyinkHzEditFieldLabel.FontSize = 15;
            app.SamplingFrequencyinkHzEditFieldLabel.FontWeight = 'bold';
            app.SamplingFrequencyinkHzEditFieldLabel.FontAngle = 'italic';
            app.SamplingFrequencyinkHzEditFieldLabel.FontColor = [0.0863 0.7216 0.8784];
            app.SamplingFrequencyinkHzEditFieldLabel.Layout.Row = 3;
            app.SamplingFrequencyinkHzEditFieldLabel.Layout.Column = 1;
            app.SamplingFrequencyinkHzEditFieldLabel.Text = 'Sampling Frequency in kHz';

            % Create SamplingFrequencyinkHzEditField
            app.SamplingFrequencyinkHzEditField = uieditfield(app.GridLayout2, 'numeric');
            app.SamplingFrequencyinkHzEditField.Limits = [0 Inf];
            app.SamplingFrequencyinkHzEditField.HorizontalAlignment = 'center';
            app.SamplingFrequencyinkHzEditField.FontName = 'Arial';
            app.SamplingFrequencyinkHzEditField.FontSize = 25;
            app.SamplingFrequencyinkHzEditField.FontAngle = 'italic';
            app.SamplingFrequencyinkHzEditField.FontColor = [1 0.4118 0.1608];
            app.SamplingFrequencyinkHzEditField.Layout.Row = 3;
            app.SamplingFrequencyinkHzEditField.Layout.Column = 2;
            app.SamplingFrequencyinkHzEditField.Value = 5;

            % Create TimeDurationofeachbitEditFieldLabel
            app.TimeDurationofeachbitEditFieldLabel = uilabel(app.GridLayout2);
            app.TimeDurationofeachbitEditFieldLabel.HorizontalAlignment = 'center';
            app.TimeDurationofeachbitEditFieldLabel.WordWrap = 'on';
            app.TimeDurationofeachbitEditFieldLabel.FontName = 'Calibri';
            app.TimeDurationofeachbitEditFieldLabel.FontSize = 15;
            app.TimeDurationofeachbitEditFieldLabel.FontWeight = 'bold';
            app.TimeDurationofeachbitEditFieldLabel.FontAngle = 'italic';
            app.TimeDurationofeachbitEditFieldLabel.FontColor = [0.0863 0.7216 0.8784];
            app.TimeDurationofeachbitEditFieldLabel.Layout.Row = 4;
            app.TimeDurationofeachbitEditFieldLabel.Layout.Column = 1;
            app.TimeDurationofeachbitEditFieldLabel.Text = 'Time Duration of each bit';

            % Create TimeDurationofeachbitEditField
            app.TimeDurationofeachbitEditField = uieditfield(app.GridLayout2, 'numeric');
            app.TimeDurationofeachbitEditField.Limits = [0 Inf];
            app.TimeDurationofeachbitEditField.ValueDisplayFormat = '%5.2f';
            app.TimeDurationofeachbitEditField.HorizontalAlignment = 'center';
            app.TimeDurationofeachbitEditField.FontName = 'Arial';
            app.TimeDurationofeachbitEditField.FontSize = 25;
            app.TimeDurationofeachbitEditField.FontAngle = 'italic';
            app.TimeDurationofeachbitEditField.FontColor = [1 0.4118 0.1608];
            app.TimeDurationofeachbitEditField.Layout.Row = 4;
            app.TimeDurationofeachbitEditField.Layout.Column = 2;
            app.TimeDurationofeachbitEditField.Value = 0.1;

            % Create SamplePointsEditFieldLabel
            app.SamplePointsEditFieldLabel = uilabel(app.GridLayout2);
            app.SamplePointsEditFieldLabel.HorizontalAlignment = 'center';
            app.SamplePointsEditFieldLabel.WordWrap = 'on';
            app.SamplePointsEditFieldLabel.FontName = 'Calibri';
            app.SamplePointsEditFieldLabel.FontSize = 15;
            app.SamplePointsEditFieldLabel.FontWeight = 'bold';
            app.SamplePointsEditFieldLabel.FontAngle = 'italic';
            app.SamplePointsEditFieldLabel.FontColor = [0.0863 0.7216 0.8784];
            app.SamplePointsEditFieldLabel.Layout.Row = 5;
            app.SamplePointsEditFieldLabel.Layout.Column = 1;
            app.SamplePointsEditFieldLabel.Text = 'Sample Points';

            % Create SamplePointsEditField
            app.SamplePointsEditField = uieditfield(app.GridLayout2, 'numeric');
            app.SamplePointsEditField.Limits = [0 Inf];
            app.SamplePointsEditField.RoundFractionalValues = 'on';
            app.SamplePointsEditField.HorizontalAlignment = 'center';
            app.SamplePointsEditField.FontName = 'Arial';
            app.SamplePointsEditField.FontSize = 25;
            app.SamplePointsEditField.FontAngle = 'italic';
            app.SamplePointsEditField.FontColor = [1 0.4118 0.1608];
            app.SamplePointsEditField.Layout.Row = 5;
            app.SamplePointsEditField.Layout.Column = 2;
            app.SamplePointsEditField.Value = 512;

            % Create NEXTButton
            app.NEXTButton = uibutton(app.GridLayout2, 'push');
            app.NEXTButton.ButtonPushedFcn = createCallbackFcn(app, @NEXTButtonPushed, true);
            app.NEXTButton.FontName = 'Arial';
            app.NEXTButton.FontSize = 30;
            app.NEXTButton.FontWeight = 'bold';
            app.NEXTButton.FontColor = [0 0 1];
            app.NEXTButton.Layout.Row = 8;
            app.NEXTButton.Layout.Column = [1 2];
            app.NEXTButton.Text = 'NEXT';

            % Create TextArea_2
            app.TextArea_2 = uitextarea(app.GridLayout2);
            app.TextArea_2.FontName = 'Times New Roman';
            app.TextArea_2.FontSize = 14;
            app.TextArea_2.FontAngle = 'italic';
            app.TextArea_2.Layout.Row = [1 2];
            app.TextArea_2.Layout.Column = [1 2];
            app.TextArea_2.Value = {'Instructions:'; '1.You can change the Parameters only before decoding. '; '2. Press NEXT to lock the values and move to Keypad.'; '3. For any information regarding any key right click on the key.'; '4. Click Encoded Graphs tab to se encoded Graphs.'; '5. Decoded Graphs will pop-up in the figure.'};

            % Create NoiseAmplitudeEditFieldLabel
            app.NoiseAmplitudeEditFieldLabel = uilabel(app.GridLayout2);
            app.NoiseAmplitudeEditFieldLabel.HorizontalAlignment = 'center';
            app.NoiseAmplitudeEditFieldLabel.WordWrap = 'on';
            app.NoiseAmplitudeEditFieldLabel.FontName = 'Calibri';
            app.NoiseAmplitudeEditFieldLabel.FontSize = 15;
            app.NoiseAmplitudeEditFieldLabel.FontWeight = 'bold';
            app.NoiseAmplitudeEditFieldLabel.FontAngle = 'italic';
            app.NoiseAmplitudeEditFieldLabel.FontColor = [0.0863 0.7216 0.8784];
            app.NoiseAmplitudeEditFieldLabel.Layout.Row = 6;
            app.NoiseAmplitudeEditFieldLabel.Layout.Column = 1;
            app.NoiseAmplitudeEditFieldLabel.Text = 'Noise Amplitude';

            % Create NoiseAmplitudeEditField
            app.NoiseAmplitudeEditField = uieditfield(app.GridLayout2, 'numeric');
            app.NoiseAmplitudeEditField.Limits = [0 Inf];
            app.NoiseAmplitudeEditField.HorizontalAlignment = 'center';
            app.NoiseAmplitudeEditField.FontName = 'Arial';
            app.NoiseAmplitudeEditField.FontSize = 25;
            app.NoiseAmplitudeEditField.FontAngle = 'italic';
            app.NoiseAmplitudeEditField.FontColor = [1 0.4118 0.1608];
            app.NoiseAmplitudeEditField.Layout.Row = 6;
            app.NoiseAmplitudeEditField.Layout.Column = 2;
            app.NoiseAmplitudeEditField.Value = 0.5;

            % Create FilterLengthEditFieldLabel
            app.FilterLengthEditFieldLabel = uilabel(app.GridLayout2);
            app.FilterLengthEditFieldLabel.HorizontalAlignment = 'center';
            app.FilterLengthEditFieldLabel.WordWrap = 'on';
            app.FilterLengthEditFieldLabel.FontName = 'Calibri';
            app.FilterLengthEditFieldLabel.FontSize = 15;
            app.FilterLengthEditFieldLabel.FontWeight = 'bold';
            app.FilterLengthEditFieldLabel.FontAngle = 'italic';
            app.FilterLengthEditFieldLabel.FontColor = [0.0863 0.7216 0.8784];
            app.FilterLengthEditFieldLabel.Layout.Row = 7;
            app.FilterLengthEditFieldLabel.Layout.Column = 1;
            app.FilterLengthEditFieldLabel.Text = 'Filter Length';

            % Create FilterLengthEditField
            app.FilterLengthEditField = uieditfield(app.GridLayout2, 'numeric');
            app.FilterLengthEditField.Limits = [0 Inf];
            app.FilterLengthEditField.RoundFractionalValues = 'on';
            app.FilterLengthEditField.HorizontalAlignment = 'center';
            app.FilterLengthEditField.FontName = 'Arial';
            app.FilterLengthEditField.FontSize = 25;
            app.FilterLengthEditField.FontAngle = 'italic';
            app.FilterLengthEditField.FontColor = [1 0.4118 0.1608];
            app.FilterLengthEditField.Layout.Row = 7;
            app.FilterLengthEditField.Layout.Column = 2;
            app.FilterLengthEditField.Value = 40;

            % Create KeyPad_Tab
            app.KeyPad_Tab = uitab(app.TabGroup);
            app.KeyPad_Tab.Title = 'KeyPad_Tab';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.KeyPad_Tab);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x'};

            % Create Button_2
            app.Button_2 = uibutton(app.GridLayout, 'push');
            app.Button_2.ButtonPushedFcn = createCallbackFcn(app, @Button_2Pushed, true);
            app.Button_2.IconAlignment = 'center';
            app.Button_2.FontName = 'Arial';
            app.Button_2.FontSize = 25;
            app.Button_2.FontWeight = 'bold';
            app.Button_2.Layout.Row = 3;
            app.Button_2.Layout.Column = 3;
            app.Button_2.Text = '2';

            % Create Button_3
            app.Button_3 = uibutton(app.GridLayout, 'push');
            app.Button_3.ButtonPushedFcn = createCallbackFcn(app, @Button_3Pushed, true);
            app.Button_3.FontName = 'Arial';
            app.Button_3.FontSize = 25;
            app.Button_3.FontWeight = 'bold';
            app.Button_3.Layout.Row = 3;
            app.Button_3.Layout.Column = 4;
            app.Button_3.Text = '3';

            % Create TextArea
            app.TextArea = uitextarea(app.GridLayout);
            app.TextArea.Editable = 'off';
            app.TextArea.HorizontalAlignment = 'right';
            app.TextArea.FontName = 'Calibri';
            app.TextArea.FontSize = 20;
            app.TextArea.FontWeight = 'bold';
            app.TextArea.Layout.Row = 2;
            app.TextArea.Layout.Column = [2 6];

            % Create Button_4
            app.Button_4 = uibutton(app.GridLayout, 'push');
            app.Button_4.ButtonPushedFcn = createCallbackFcn(app, @Button_4Pushed, true);
            app.Button_4.FontName = 'Arial';
            app.Button_4.FontSize = 25;
            app.Button_4.FontWeight = 'bold';
            app.Button_4.Layout.Row = 4;
            app.Button_4.Layout.Column = 2;
            app.Button_4.Text = '4';

            % Create Button_5
            app.Button_5 = uibutton(app.GridLayout, 'push');
            app.Button_5.ButtonPushedFcn = createCallbackFcn(app, @Button_5Pushed, true);
            app.Button_5.FontName = 'Arial';
            app.Button_5.FontSize = 25;
            app.Button_5.FontWeight = 'bold';
            app.Button_5.Layout.Row = 4;
            app.Button_5.Layout.Column = 3;
            app.Button_5.Text = '5';

            % Create Button_6
            app.Button_6 = uibutton(app.GridLayout, 'push');
            app.Button_6.ButtonPushedFcn = createCallbackFcn(app, @Button_6Pushed, true);
            app.Button_6.FontName = 'Arial';
            app.Button_6.FontSize = 25;
            app.Button_6.FontWeight = 'bold';
            app.Button_6.Layout.Row = 4;
            app.Button_6.Layout.Column = 4;
            app.Button_6.Text = '6';

            % Create Button_7
            app.Button_7 = uibutton(app.GridLayout, 'push');
            app.Button_7.ButtonPushedFcn = createCallbackFcn(app, @Button_7Pushed, true);
            app.Button_7.FontName = 'Arial';
            app.Button_7.FontSize = 25;
            app.Button_7.FontWeight = 'bold';
            app.Button_7.Layout.Row = 5;
            app.Button_7.Layout.Column = 2;
            app.Button_7.Text = '7';

            % Create Button_8
            app.Button_8 = uibutton(app.GridLayout, 'push');
            app.Button_8.ButtonPushedFcn = createCallbackFcn(app, @Button_8Pushed, true);
            app.Button_8.FontName = 'Arial';
            app.Button_8.FontSize = 25;
            app.Button_8.FontWeight = 'bold';
            app.Button_8.Layout.Row = 5;
            app.Button_8.Layout.Column = 3;
            app.Button_8.Text = '8';

            % Create Button_9
            app.Button_9 = uibutton(app.GridLayout, 'push');
            app.Button_9.ButtonPushedFcn = createCallbackFcn(app, @Button_9Pushed, true);
            app.Button_9.FontName = 'Arial';
            app.Button_9.FontSize = 25;
            app.Button_9.FontWeight = 'bold';
            app.Button_9.Layout.Row = 5;
            app.Button_9.Layout.Column = 4;
            app.Button_9.Text = '9';

            % Create Button_10
            app.Button_10 = uibutton(app.GridLayout, 'push');
            app.Button_10.ButtonPushedFcn = createCallbackFcn(app, @Button_10Pushed, true);
            app.Button_10.VerticalAlignment = 'bottom';
            app.Button_10.FontName = 'Arial';
            app.Button_10.FontSize = 40;
            app.Button_10.Layout.Row = 6;
            app.Button_10.Layout.Column = 2;
            app.Button_10.Text = '*';

            % Create Button_11
            app.Button_11 = uibutton(app.GridLayout, 'push');
            app.Button_11.ButtonPushedFcn = createCallbackFcn(app, @Button_11Pushed, true);
            app.Button_11.FontName = 'Arial';
            app.Button_11.FontSize = 25;
            app.Button_11.FontWeight = 'bold';
            app.Button_11.Layout.Row = 6;
            app.Button_11.Layout.Column = 3;
            app.Button_11.Text = '0';

            % Create Button_12
            app.Button_12 = uibutton(app.GridLayout, 'push');
            app.Button_12.ButtonPushedFcn = createCallbackFcn(app, @Button_12Pushed, true);
            app.Button_12.FontName = 'Arial';
            app.Button_12.FontSize = 25;
            app.Button_12.FontWeight = 'bold';
            app.Button_12.Layout.Row = 6;
            app.Button_12.Layout.Column = 4;
            app.Button_12.Text = '#';

            % Create KEYPADLabel
            app.KEYPADLabel = uilabel(app.GridLayout);
            app.KEYPADLabel.HorizontalAlignment = 'center';
            app.KEYPADLabel.FontName = 'Arial';
            app.KEYPADLabel.FontSize = 30;
            app.KEYPADLabel.FontWeight = 'bold';
            app.KEYPADLabel.FontColor = [0.851 0.3255 0.098];
            app.KEYPADLabel.Layout.Row = 1;
            app.KEYPADLabel.Layout.Column = [2 6];
            app.KEYPADLabel.Text = 'KEYPAD';

            % Create ClearButton
            app.ClearButton = uibutton(app.GridLayout, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.WordWrap = 'on';
            app.ClearButton.FontName = 'Arial';
            app.ClearButton.FontSize = 18;
            app.ClearButton.FontWeight = 'bold';
            app.ClearButton.Layout.Row = 3;
            app.ClearButton.Layout.Column = 6;
            app.ClearButton.Text = 'Clear';

            % Create EnterButton
            app.EnterButton = uibutton(app.GridLayout, 'push');
            app.EnterButton.ButtonPushedFcn = createCallbackFcn(app, @EnterButtonPushed, true);
            app.EnterButton.WordWrap = 'on';
            app.EnterButton.FontName = 'Arial';
            app.EnterButton.FontSize = 18;
            app.EnterButton.FontWeight = 'bold';
            app.EnterButton.Layout.Row = [5 6];
            app.EnterButton.Layout.Column = 6;
            app.EnterButton.Text = 'Enter';

            % Create AButton
            app.AButton = uibutton(app.GridLayout, 'push');
            app.AButton.ButtonPushedFcn = createCallbackFcn(app, @AButtonPushed, true);
            app.AButton.FontName = 'Arial';
            app.AButton.FontSize = 25;
            app.AButton.FontWeight = 'bold';
            app.AButton.Layout.Row = 3;
            app.AButton.Layout.Column = 5;
            app.AButton.Text = 'A';

            % Create BButton
            app.BButton = uibutton(app.GridLayout, 'push');
            app.BButton.ButtonPushedFcn = createCallbackFcn(app, @BButtonPushed, true);
            app.BButton.FontName = 'Arial';
            app.BButton.FontSize = 25;
            app.BButton.FontWeight = 'bold';
            app.BButton.Layout.Row = 4;
            app.BButton.Layout.Column = 5;
            app.BButton.Text = 'B';

            % Create CButton
            app.CButton = uibutton(app.GridLayout, 'push');
            app.CButton.ButtonPushedFcn = createCallbackFcn(app, @CButtonPushed, true);
            app.CButton.FontName = 'Arial';
            app.CButton.FontSize = 25;
            app.CButton.FontWeight = 'bold';
            app.CButton.Layout.Row = 5;
            app.CButton.Layout.Column = 5;
            app.CButton.Text = 'C';

            % Create DButton
            app.DButton = uibutton(app.GridLayout, 'push');
            app.DButton.ButtonPushedFcn = createCallbackFcn(app, @DButtonPushed, true);
            app.DButton.FontName = 'Arial';
            app.DButton.FontSize = 25;
            app.DButton.FontWeight = 'bold';
            app.DButton.Layout.Row = 6;
            app.DButton.Layout.Column = 5;
            app.DButton.Text = 'D';

            % Create BackButton
            app.BackButton = uibutton(app.GridLayout, 'push');
            app.BackButton.ButtonPushedFcn = createCallbackFcn(app, @BackButtonPushed, true);
            app.BackButton.WordWrap = 'on';
            app.BackButton.FontName = 'Arial';
            app.BackButton.FontSize = 18;
            app.BackButton.FontWeight = 'bold';
            app.BackButton.Layout.Row = 4;
            app.BackButton.Layout.Column = 6;
            app.BackButton.Text = 'Back';

            % Create Button_1
            app.Button_1 = uibutton(app.GridLayout, 'push');
            app.Button_1.ButtonPushedFcn = createCallbackFcn(app, @Button_1Pushed, true);
            app.Button_1.IconAlignment = 'center';
            app.Button_1.FontName = 'Arial';
            app.Button_1.FontSize = 25;
            app.Button_1.FontWeight = 'bold';
            app.Button_1.Layout.Row = 3;
            app.Button_1.Layout.Column = 2;
            app.Button_1.Text = '1';

            % Create EncodedGraphsTab
            app.EncodedGraphsTab = uitab(app.TabGroup);
            app.EncodedGraphsTab.Title = 'Encoded Graphs';

            % Create UIAxes
            app.UIAxes = uiaxes(app.EncodedGraphsTab);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [0 318 561 212];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.EncodedGraphsTab);
            title(app.UIAxes2, 'Title')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Position = [0 32 560 201];

            % Create Result
            app.Result = uitextarea(app.UIFigure);
            app.Result.Editable = 'off';
            app.Result.HorizontalAlignment = 'center';
            app.Result.FontName = 'Times New Roman';
            app.Result.FontSize = 22;
            app.Result.FontWeight = 'bold';
            app.Result.FontColor = [0.4941 0.1843 0.5569];
            app.Result.Position = [144 571 560 68];

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create ClearstheentirescreenanderasestheprocessMenu
            app.ClearstheentirescreenanderasestheprocessMenu = uimenu(app.ContextMenu);
            app.ClearstheentirescreenanderasestheprocessMenu.Text = 'Clears the entire screen and erases the process';
            
            % Assign app.ContextMenu
            app.ClearButton.ContextMenu = app.ContextMenu;

            % Create ContextMenu2
            app.ContextMenu2 = uicontextmenu(app.UIFigure);

            % Create BackspacedeletesthedatacorrespondingtolastinputkeyMenu
            app.BackspacedeletesthedatacorrespondingtolastinputkeyMenu = uimenu(app.ContextMenu2);
            app.BackspacedeletesthedatacorrespondingtolastinputkeyMenu.Text = 'Backspace, deletes the data corresponding to last input key';
            
            % Assign app.ContextMenu2
            app.BackButton.ContextMenu = app.ContextMenu2;

            % Create ContextMenu3
            app.ContextMenu3 = uicontextmenu(app.UIFigure);

            % Create Menu
            app.Menu = uimenu(app.ContextMenu3);
            app.Menu.Text = 'Press to make the sequence entered as the encoded sequence and move to decoding';
            
            % Assign app.ContextMenu3
            app.EnterButton.ContextMenu = app.ContextMenu3;

            % Create ContextMenu4
            app.ContextMenu4 = uicontextmenu(app.UIFigure);

            % Create FixtheabovevaluesandstartdecodingMenu
            app.FixtheabovevaluesandstartdecodingMenu = uimenu(app.ContextMenu4);
            app.FixtheabovevaluesandstartdecodingMenu.Text = 'Fix the above values and start decoding';
            
            % Assign app.ContextMenu4
            app.NEXTButton.ContextMenu = app.ContextMenu4;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Expt3_Keypad_app_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end