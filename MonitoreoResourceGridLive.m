function MonitoreoVivoResourceGrid
    % Configuración inicial y parámetros de SDR/banda
    clear; clc;
    radioOptions = hSDRBase.getDeviceNameOptions;
    rx = hSDRReceiver(radioOptions{10});
    antennaOptions = getAntennaOptions(rx);
    rx.ChannelMapping = antennaOptions(1);
    rx.Gain = 50;
    band = "n78";
    GSCN = 7929;
    rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN);
    scs = "30kHz";
    nrbSSB = 30;
    scsNumeric = double(extract(scs,digitsPattern));
    ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric);
    rx.SampleRate = ofdmInfo.SampleRate;
    
    % Parámetros de captura y visualización
    captureDuration = seconds(0.03);  % Duración de cada captura (~3 frames de 10ms cada uno)
    pauseInterval = 0.1;              % Pausa entre iteraciones en segundos
    monitoreoTiempo = 600*2;           % Duración total monitorización en segundos (e.g. 20 min)
    
    % Inicializar la figura para mostrar el resource grid
    hFig = figure('Name','Monitoreo Vivo Resource Grid 5G','NumberTitle','off');
    hImg = imagesc(zeros(360,56)); % Ajustar tamaño según configuración (filas=subcarriers, cols=símbolos OFDM)
    axis xy; colormap('jet'); colorbar;
    xlabel('Símbolos OFDM'); ylabel('Subcarriers');
    title('Resource Grid en vivo 5G');
    hold on;
    hold off;
    
    % Iniciar temporizador
    tic;
    
    while ishandle(hFig) && toc < monitoreoTiempo
        try
            waveform = capture(rx, captureDuration);
            % Procesar la señal con función findSSB adaptada
            [detectedSSB, gridSSB, burstRect, burstText, ncellid] = findSSBrobusto(waveform,...
                rx.CenterFrequency, scs, rx.SampleRate);
            
            % Solo actualizar si grid válido y detectado
            if detectedSSB && ~isempty(gridSSB)
                set(hImg, 'CData', abs(gridSSB));
                
          
                title(sprintf('Resource Grid t=%.2fs | CellID=%d', toc, ncellid));
                drawnow;
            else
                % Opcional: mostrar mensaje si no se detecta SSB
                title(sprintf('No se detectó SSB en t=%.2fs', toc));
                drawnow;
            end
        catch ME
            warning('Error en procesamiento: %s', ME.message);
            % Continuar sin detener la captura
        end
        pause(pauseInterval); % Controla la tasa de actualización
    end
    
    % Limpieza final
    if ishandle(hFig), close(hFig); end
    release(rx);
    fprintf('Monitoreo finalizado.\n');
end

function [detectedSSB, gridSSB, rectSSB, burstText, ncellid] = findSSBrobusto(waveform, centerFrequency, scs, sampleRate)
    % Parámetros de ayuda
    try
        ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs, centerFrequency);
        scsNumeric = double(extract(scs, digitsPattern));
        searchBW = 3*scsNumeric;
        displayFigure = false;
        
        % Corrección de frecuencia
        [correctedWaveform,~,NID2] = hSSBurstFrequencyCorrect(waveform, ssbBlockPattern, sampleRate, searchBW, displayFigure);
        
        nrbSSB = 20;
        refGrid = zeros([nrbSSB*12 2]);
        refGrid(nrPSSIndices,2) = nrPSS(NID2);
        nSlot = 0;
        
        % Estimación de timing
        timingOffset = nrTimingEstimate(correctedWaveform,nrbSSB,scsNumeric,nSlot,refGrid,SampleRate=sampleRate);
        correctedWaveform = correctedWaveform(1+timingOffset:end,:);
        
        % Demodulación OFDM
        rxGrid = nrOFDMDemodulate(correctedWaveform,nrbSSB,scsNumeric,nSlot,SampleRate=sampleRate);
        
        % Validar tamaño para acceder a símbolo 2:5
        numSymbols = size(rxGrid, 2);
        if numSymbols < 5
            detectedSSB = false;
            gridSSB = [];
            rectSSB = [];
            burstText = '';
            ncellid = -1;
            return;
        end
        
        rxGridSSB = rxGrid(:, 2:5, :); % Símbolos OFDM centrales
        
        % Correlación estimación NID1 usando SSS
        sssIndices = nrSSSIndices;
        sssRx = nrExtractResources(sssIndices, rxGridSSB);
        
        sssEst = zeros(1,336);
        for NID1 = 0:335
            ncellidTemp = (3*NID1) + NID2;
            sssRef = nrSSS(ncellidTemp);
            sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef),1)).^2);
        end
        
        [maxCorr, idx] = max(sssEst);
        NID1 = idx - 1;
        ncellid = (3*NID1) + NID2;
        
        % Estimar burst SSB más fuerte usando DM-RS
        dmrsIndices = nrPBCHDMRSIndices(ncellid);
        dmrsEst = zeros(1,8);
        for ibar_SSB = 0:7
            refGridDMRS = zeros([240 4]);
            refGridDMRS(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
            [hest,nest] = nrChannelEstimate(rxGridSSB,refGridDMRS,'AveragingWindow',[0 1]);
            dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);
        end
        
        [~, ibar_SSB] = max(dmrsEst);
        strongestSSBIdx = ibar_SSB - 1;
        detectedSSB = maxCorr > 1e-3;
        burstText = sprintf('Strongest SSB: %d', strongestSSBIdx);
        
        % Construir el resource grid para despliegue
        demodRB = 30;
        gridSSB = nrOFDMDemodulate(correctedWaveform, demodRB, scsNumeric, nSlot, SampleRate=sampleRate);
        last = min(56, size(gridSSB, 2));
        gridSSB = gridSSB(:, 1:last, 1);
        
        ssbFreqOrigin = 12*(demodRB-nrbSSB)/2 + 1;
        startSymbol = 1;
        numSymbolsSSB = 4;
        rectSSB = [startSymbol+0.5, ssbFreqOrigin-0.5, numSymbolsSSB, 12*nrbSSB];
        
    catch ME_inner
        % Fallo procesamiento: devolver vacíos y marcado como no detectado
        detectedSSB = false;
        gridSSB = [];
        rectSSB = [];
        burstText = sprintf('Error: %s', ME_inner.message);
        ncellid = -1;
    end
end
