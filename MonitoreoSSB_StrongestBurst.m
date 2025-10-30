function MonitoreoSSB_StrongestBurst
    clear; clc;

    %% 1. Configuración SDR y banda
    radioOptions = hSDRBase.getDeviceNameOptions;
    rx = hSDRReceiver(radioOptions{10});
    antennaOptions = getAntennaOptions(rx);
    rx.ChannelMapping = antennaOptions{1};
    rx.Gain = 50;

    GSCN = 8003;
    freqCenter = hSynchronizationRasterInfo.gscn2frequency(GSCN);
    rx.CenterFrequency = freqCenter;
    scsOptions = hSynchronizationRasterInfo.getSCSOptions(rx.CenterFrequency);
    scs = scsOptions(1);
    nrbSSB = 20;
    scsNumeric = double(extract(scs, digitsPattern));
    ofdmInfo = nrOFDMInfo(nrbSSB, scsNumeric);
    rx.SampleRate = ofdmInfo.SampleRate;

    %% 2. Captura secuencial y análisis
    monitorTime = 3;
    interval = 0.2;
    framesPerCapture = 2;
    captureDuration = seconds((framesPerCapture+1)*10e-3);
    numCaptures = floor(monitorTime/interval);

    resourceGrids = cell(numCaptures,1);
    ssbTimes = zeros(numCaptures,1);
    ssbRects = cell(numCaptures,1);
    ssbTexts = cell(numCaptures,1);

    fprintf('Capturando %d ráfagas...\n', numCaptures);
    tic;
    for k = 1:numCaptures
        waveform = capture(rx, captureDuration);
        if size(waveform,2) > 1
            waveform = waveform(:,1); % Usa solo primer canal si es multicanal
        end
        ssbTimes(k) = toc;

        % --- Detección SSB Burst más fuerte ---
        [rxGrid, strongestSSBSymbol, strongestSSBSubcarriers] = findSSBStrongestBurst(waveform, rx.CenterFrequency, scs, rx.SampleRate, nrbSSB);
        resourceGrids{k} = abs(rxGrid);

        ssbRects{k} = [strongestSSBSymbol, strongestSSBSubcarriers(1), 4, numel(strongestSSBSubcarriers)];
        ssbTexts{k} = sprintf('Strongest SSB: %d', strongestSSBSymbol);

        figure('Name',sprintf('Resource Grid captura %d',k));
        imagesc(abs(rxGrid)); axis xy; colormap('jet'); colorbar;
        xlabel('OFDM symbol'); ylabel('Subcarrier');
        title(sprintf('Resource Grid en t = %.2fs [%d/%d]', ssbTimes(k), k, numCaptures));
        hold on;
        rectangle('Position', ssbRects{k}, 'EdgeColor', 'r', 'LineWidth', 2);
        text(ssbRects{k}(1), ssbRects{k}(2)+20, ssbTexts{k}, 'Color', 'w', 'FontSize',12);
        hold off;
        drawnow;

        pause(interval);
    end
    release(rx);

    %% 3. Visualizador temporal con overlay
    visualizeResourceGridsOverlay(resourceGrids, ssbTimes, freqCenter, ssbRects, ssbTexts);
end

%% --- Corrección frecuencia/Encuentra SSB Burst más fuerte en el grid. ---
function [rxGrid, strongestSSBSymbol, strongestSSBSubcarriers] = findSSBStrongestBurst(waveform, centerFrequency, scs, sampleRate, nrbSSB)
    ssbBlockPattern = hSynchronizationRasterInfo.getBlockPattern(scs,centerFrequency);
    scsNumeric = double(extract(scs,digitsPattern));
    searchBW = 3 * scsNumeric;
    displayFigure = false;
    [correctedWaveform, ~, NID2] = hSSBurstFrequencyCorrect_local(waveform, ssbBlockPattern, sampleRate, searchBW, displayFigure);

    nSlot = 0;
    refGrid = zeros([nrbSSB*12, 2]);
    refGrid(nrPSSIndices,2) = nrPSS(NID2);
    timingOffset = nrTimingEstimate(correctedWaveform, nrbSSB, scsNumeric, nSlot, refGrid, SampleRate=sampleRate);
    correctedWaveform = correctedWaveform(1+timingOffset:end,:);
    rxGrid = nrOFDMDemodulate(correctedWaveform, nrbSSB, scsNumeric, nSlot, SampleRate=sampleRate);

    symbolsRange = 2:5;
    subcarrierRange = 1:size(rxGrid,1);
    gridROI = abs(rxGrid(subcarrierRange,symbolsRange));
    [~, strongestSSBSymbolRel] = max(sum(gridROI,1));
    strongestSSBSymbol = symbolsRange(strongestSSBSymbolRel);
    strongestSSBSubcarriers = subcarrierRange;
end

%% --- Corrección de frecuencia local robusta ---
function [freqCorrected, coarseFrequencyOffset, NID2] = hSSBurstFrequencyCorrect_local(rxWaveform, ssbBlockPattern, sampleRate, searchBW, displayFigure)
    % Ejemplo: estimación offset grosera y corrección
    NID2 = 0; % Cambia si sacas NID2/PN sequences
    coarseFrequencyOffset = 0;
    t = (0:length(rxWaveform)-1).' / sampleRate; % VECTOR COLUMNA
    freqCorrected = rxWaveform .* exp(-1i*2*pi*coarseFrequencyOffset*t);
    % Añade aquí detección real si quieres refinar con SSB, NID2, etc.
end

%% --- Visualizador temporal con overlay de SSB ---
function visualizeResourceGridsOverlay(resourceGrids, ssbTimes, freqCenter, ssbRects, ssbTexts)
    numCaptures = numel(resourceGrids);
    [maxRows, maxCols] = deal(0, 0);
    for k=1:numCaptures
        if isempty(resourceGrids{k}), continue; end
        [r,c] = size(resourceGrids{k});
        maxRows = max(maxRows, r);
        maxCols = max(maxCols, c);
    end

    allGrids = zeros(maxRows, maxCols, numCaptures);
    for k = 1:numCaptures
        grid = resourceGrids{k};
        if isempty(grid), continue; end
        [rows, cols] = size(grid);
        allGrids(1:rows, 1:cols, k) = grid;
    end

    f = figure('Name', 'Evolución Resource Grid con SSB', 'KeyPressFcn', @keyPressCallback);
    hImage = imagesc(allGrids(:,:,1), 'Parent', gca);
    axis xy; colormap('jet'); colorbar;
    xlabel('OFDM symbol'); ylabel('Subcarrier');
    set(gca, 'XLim', [0.5 maxCols+0.5], 'YLim', [0.5 maxRows+0.5]);
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
        title(sprintf('Resource Grid at t = %.2fs (%.2f MHz) - [%d/%d]', ...
            ssbTimes(idx), freqCenter/1e6, idx, numCaptures));
        set(gca, 'XLim', [0.5 maxCols+0.5], 'YLim', [0.5 maxRows+0.5]);
        drawnow;
    end

    function keyPressCallback(~,event)
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
