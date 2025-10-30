function [waveformsAll, resourceGrids, ssbTimes, powerVec, snrVec, cellIDVec] = MonitoreoContinuo2
    clear; clc;
    % --- Configuración SDR, banda y parámetros ---
    radioOptions = hSDRBase.getDeviceNameOptions;
    rx = hSDRReceiver(radioOptions{10}); % o B210 si tu SDR tiene ese nombre en radioOptions
    antennaOptions = getAntennaOptions(rx);
    rx.ChannelMapping = antennaOptions(1);
    rx.Gain = 50;
    band = "n78";
    GSCN = 8003;
    rx.CenterFrequency = hSynchronizationRasterInfo.gscn2frequency(GSCN);
    scs = "30kHz";
    nrbSSB = 20;
    scsNumeric = double(extract(scs, digitsPattern));
    ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric);
    rx.SampleRate = ofdmInfo.SampleRate;
    % --- Parámetros de captura ---
    monitorTime = 1; % segundos
    interval = 0.1;  % segundos
    framesPerCapture = 2;
    captureDuration = seconds((framesPerCapture+1)*10e-3);
    numCaptures = floor(monitorTime/interval);
    % --- Captura secuencial ---
    resourceGrids = cell(numCaptures,1);
    ssbTimes = zeros(numCaptures,1);
    powerVec = zeros(numCaptures,1);
    snrVec = zeros(numCaptures,1);
    cellIDVec = zeros(numCaptures,1);
    ssbRects = cell(numCaptures,1);
    ssbTexts = cell(numCaptures,1);
    waveformsAll = cell(numCaptures,1); % Aquí guardamos todas las formas de onda
    fprintf('Capturando %d ráfagas...\n', numCaptures);
    tic;
    for k = 1:numCaptures
        waveform = capture(rx, captureDuration);
        waveformsAll{k} = waveform; % Guardar la forma de onda
        ssbTimes(k) = toc;
        % --- Detección SSB usando findSSB (como en el ejemplo) ---
        [detectedSSB, gridSSB, burstRect, burstText, ncellid] = findSSB(waveform, rx.CenterFrequency, scs, rx.SampleRate);
        powerVec(k) = 10*log10(mean(abs(waveform).^2));
        snrVec(k) = estimateSNR(waveform);
        cellIDVec(k) = ncellid;
        resourceGrids{k} = gridSSB;
        ssbRects{k} = burstRect;
        ssbTexts{k} = burstText;
        fprintf('[%.2fs] Potencia=%.1f dB | SNR=%.1f dB | cellID=%d\n', ssbTimes(k), powerVec(k), snrVec(k), ncellid);
        % Mostrar figura individual estilo ejemplo
              % Mostrar figura individual estilo ejemplo (comentado)
        % figure('Name',sprintf('Resource Grid captura %d',k));
        % imagesc(abs(gridSSB)); axis xy; colormap('jet'); colorbar;
        % xlabel('OFDM symbol'); ylabel('Subcarrier');
        % title(sprintf('Resource Grid en t = %.2fs [%d/%d]', ssbTimes(k), k, numCaptures));
        % if detectedSSB
        %     hold on;
        %     rectangle('Position', burstRect, 'EdgeColor', 'r', 'LineWidth', 2);
        %     text(burstRect(1), burstRect(2)+20, burstText, 'Color', 'w', 'FontSize',12);
        %     hold off;
        % end
        % drawnow;
        pause(interval);
    end
    release(rx);
    % --- Visualizador interactivo temporal (slider) con overlay ---
    visualizeResourceGridsOverlay(resourceGrids, ssbTimes, rx.CenterFrequency, ssbRects, ssbTexts);
end


function [detectedSSB, gridSSB, rectSSB, burstText, ncellid] = findSSB(waveform,centerFrequency,scs,sampleRate)
    ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs,centerFrequency);
    scsNumeric = double(extract(scs,digitsPattern));
    searchBW = 3*scsNumeric;
    displayFigure = false;
    [correctedWaveform,~,NID2] = hSSBurstFrequencyCorrect(waveform,ssbBlockPattern,sampleRate,searchBW,displayFigure);
    nrbSSB = 20;
    refGrid = zeros([nrbSSB*12 2]);
    refGrid(nrPSSIndices,2) = nrPSS(NID2);
    nSlot = 0;
    timingOffset = nrTimingEstimate(correctedWaveform,nrbSSB,scsNumeric,nSlot,refGrid,SampleRate=sampleRate);
    correctedWaveform = correctedWaveform(1+timingOffset:end,:);
    rxGrid = nrOFDMDemodulate(correctedWaveform,nrbSSB,scsNumeric,nSlot,SampleRate=sampleRate);

    % Extrae los símbolos centrales donde está el SSB (tipicamente 2:5)
    rxGridSSB = rxGrid(:,2:5,:);
    sssIndices = nrSSSIndices;
    sssRx = nrExtractResources(sssIndices,rxGridSSB);

    % Correlación para estimar NID1 usando SSS
    sssEst = zeros(1,336);
    for NID1 = 0:335
        ncellidTemp = (3*NID1) + NID2;
        sssRef = nrSSS(ncellidTemp);
        sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef),1)).^2);
    end
    [maxCorr, idx] = max(sssEst);
    NID1 = idx - 1;
    ncellid = (3*NID1) + NID2;
    % Busca el burst SSB más fuerte
    dmrsIndices = nrPBCHDMRSIndices(ncellid);
    dmrsEst = zeros(1,8);
    for ibar_SSB = 0:7
        refGrid = zeros([240 4]);
        refGrid(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
        [hest,nest] = nrChannelEstimate(rxGridSSB,refGrid,'AveragingWindow',[0 1]);
        dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);
    end
    [~, ibar_SSB] = max(dmrsEst);
    strongestSSBIdx = ibar_SSB - 1;
    detectedSSB = maxCorr > 1e-3;
    burstText = sprintf('Strongest SSB: %d', strongestSSBIdx);

    % RECURSO GRID Y RECTÁNGULO ESTILO EJEMPLO
    demodRB = 30;
    gridSSB = nrOFDMDemodulate(correctedWaveform,demodRB,scsNumeric,nSlot,SampleRate=sampleRate);
    last = min(56,size(gridSSB,2));
    gridSSB = gridSSB(:,1:last,1);
    ssbFreqOrigin = 12*(demodRB-nrbSSB)/2 + 1;
    startSymbol = 1;
    numSymbolsSSB = 4;
    rectSSB = [startSymbol+0.5 ssbFreqOrigin-0.5 numSymbolsSSB 12*nrbSSB];
end

function visualizeResourceGridsOverlay(resourceGrids, ssbTimes, freqCenter, ssbRects, ssbTexts)
    numCaptures = numel(resourceGrids);
    fixedRows = size(resourceGrids{1},1);
    fixedCols = size(resourceGrids{1},2);
    allGrids = zeros(fixedRows, fixedCols, numCaptures);
    for k = 1:numCaptures
        grid = resourceGrids{k};
        if isempty(grid), continue; end
        allGrids(1:size(grid,1), 1:size(grid,2), k) = abs(grid);
    end
    f = figure('Name','Evolución Resource Grid con SSB','KeyPressFcn',@keyPressCallback);
    hImage = imagesc(allGrids(:,:,1), 'Parent', gca);
    axis xy; colormap('jet'); colorbar;
    xlabel('OFDM symbol'); ylabel('Subcarrier');
    set(gca, 'XLim', [0.5 fixedCols+0.5], 'YLim', [0.5 fixedRows+0.5]);
    title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [1/%d]', ssbTimes(1), freqCenter/1e6, numCaptures));
    hSlider = uicontrol('Style','slider','Min',1,'Max',numCaptures,'Value',1,...
        'SliderStep', [1/(numCaptures-1), 1/(numCaptures-1)],'Position', [150 20 300 20]);
    hSlider.Callback = @(src,~) updateGrid(round(src.Value));
    overlayRect = rectangle('Position', ssbRects{1}, 'EdgeColor', 'r', 'LineWidth', 2);
    overlayText = text(ssbRects{1}(1), ssbRects{1}(2)+20, ssbTexts{1}, 'Color','w','FontSize',12);
    function updateGrid(idx)
        set(hImage, 'CData', allGrids(:,:,idx));
        set(overlayRect, 'Position', ssbRects{idx});
        set(overlayText, 'Position', [ssbRects{idx}(1), ssbRects{idx}(2)+20]);
        set(overlayText, 'String', ssbTexts{idx});
        title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [%d/%d]', ssbTimes(idx), freqCenter/1e6, idx, numCaptures));
        set(gca, 'XLim', [0.5 fixedCols+0.5], 'YLim', [0.5 fixedRows+0.5]); 
        drawnow;
    end
    function keyPressCallback(~, event)
        idx = round(hSlider.Value);
        switch event.Key
            case 'rightarrow'
                idx = min(idx+1, numCaptures);
            case 'leftarrow'
                idx = max(idx-1, 1);
            otherwise
                return;
        end
        hSlider.Value = idx;
        updateGrid(idx);
    end
end

function SNRdB = estimateSNR(signal)
    noiseEst = median(abs(signal)).^2;
    signalEst = mean(abs(signal).^2);
    SNRdB = 10*log10(signalEst/noiseEst+eps);
end
