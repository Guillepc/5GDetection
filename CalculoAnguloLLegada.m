function [angulosRecepcion, diferenciasMediasFase, diferenciasFasePorCaptura, waveformsAll, tiemposCaptura, lambda] = MonitoreoConAngulo()
    clear; clc; 
    % Configuración SDR
    radioOptions = hSDRBase.getDeviceNameOptions; 
    rx = hSDRReceiver(radioOptions{10}); % USRP B210
    rx.ChannelMapping = [1 2]; % Usa dos antenas simultáneamente
    rx.Gain = 50;
    
    band = "n78"; 
    GSCN = 7929;
    rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN); 
    scs = "30kHz"; 
    nrbSSB = 15; 
    scsNumeric = double(extract(scs, digitsPattern)); 
    ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric); 
    rx.SampleRate = ofdmInfo.SampleRate; 
    
    monitorTime = 1; interval = 0.2; framesPerCapture = 1;
    captureDuration = seconds((framesPerCapture+1)*10e-3);
    numCaptures = floor(monitorTime/interval);
    
    angulosRecepcion = zeros(numCaptures,1);
    diferenciasMediasFase = zeros(numCaptures,1);
    diferenciasFasePorCaptura = cell(numCaptures,1);
    waveformsAll = cell(numCaptures,1);
    tiemposCaptura = zeros(numCaptures,1);
    
    distanciaAntenas = 0.016; % 1.6 cm
    c = 3e8; % velocidad de la luz
    lambda = c/rx.CenterFrequency;
    
    fprintf('Capturando %d señales con dos antenas...\n', numCaptures);
    tic;
    
    for k = 1:numCaptures
        waveform = capture(rx, captureDuration); % waveforms N x 2 (2 antenas)
        waveformsAll{k} = waveform;
        tiemposCaptura(k) = toc;
        
        IQ1 = waveform(:,1);
        IQ2 = waveform(:,2);
        
        % Extraer fase
        phi1 = angle(IQ1);
        phi2 = angle(IQ2);
        
        % Diferencia de fase para cada muestra
        deltaPhiPorMuestra = phi2 - phi1;
        diferenciasFasePorCaptura{k} = deltaPhiPorMuestra;
        
        % Diferencia media de fase
        deltaPhiMedia = mean(deltaPhiPorMuestra);
        diferenciasMediasFase(k) = deltaPhiMedia;
        
        % Calcular ángulo de llegada
        sinTheta = deltaPhiMedia*lambda/(2*pi*distanciaAntenas);
        sinTheta = min(max(sinTheta, -1), 1);
        theta = asin(sinTheta);
        
        angulosRecepcion(k) = rad2deg(theta);
        
        fprintf('[%.2fs] Ángulo de llegada = %.2f grados\n', tiemposCaptura(k), angulosRecepcion(k));
        pause(interval);
    end
    
    release(rx);
end

[angulos, diferenciasMedias, diferenciasPorCaptura, waveforms, tiempos, lambda] = MonitoreoConAngulo();
