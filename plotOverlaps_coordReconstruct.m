function plotOverlaps_coordReconstruct(initialBeamParams0,kickStrengths, offsets, planePos, planeOffsets)
    %PLOTOVERLAPS Plot "position scan" plots with variable kickStrengths
    % - initialBeamParams0: List of initial beam parameters for each kick strength
    % - offsets: offset scan array [mm]
    %
    % Note: Structure offset is treated internally in kickSim
    %
    % Assumption: Only a single structure (i.e. negative planePos)!
    %
    % Example run:
    % plotOverlaps_coordReconstruct([[0,0]', [.3,0.2]', [-.5,0.0]', [0.2,-0.5]'], [0,0.5,1.0,1.5], -2:.5:2, [1,2,-3,4],[-1,1,-1,1])
    
    if ~exist('planeOffsets', 'var')
        planeOffsets = zeros(1,length(planePos));
    end
    
    %% Sanity checks
    assert( length(initialBeamParams0) == length(kickStrengths) )
    assert( length(planePos) == length(planeOffsets) )
    
    planePos_structureIdx = -1;
    for i=1:length(planePos)
        if planePos(i) < 0
            assert (planePos_structureIdx==-1) %There can only be one!
            planePos_structureIdx = i;
        end
    end
    
    %% Side view plots
    figure(10);
    clf();
    hold on;
    colorCycle=colormap('lines');
    refPlotHandles = [];
    
    %% Simulate
    hitPositions     = zeros(length(kickStrengths),length(planePos),length(offsets));
    hitPositions_ref = zeros(length(kickStrengths),length(planePos));
    for kickStrengthIdx=1:length(kickStrengths)
        initialBeamParams = ones(2,length(offsets)) .* initialBeamParams0(:,kickStrengthIdx);
        initialBeamParams(1,:) = initialBeamParams(1,:)+offsets;
        
        kickStrength = ones(1,length(offsets)) .* kickStrengths(kickStrengthIdx);
        
        %The scan
        trueBeamParams = kickSim(initialBeamParams,kickStrength, planePos, ...
            false, planeOffsets(planePos_structureIdx) );
        %Reference beam, single beam with offset=0
        trueBeamParams_ref = kickSim(initialBeamParams0(:,kickStrengthIdx),kickStrength(kickStrengthIdx), planePos, ...
            false, planeOffsets(planePos_structureIdx));
        
        %Make shifted hit positions
        for planeIdx=1:length(planePos)
            hitPositions(kickStrengthIdx,planeIdx,:) = trueBeamParams(planeIdx+1,:,1) + planeOffsets(planeIdx);
            hitPositions_ref(kickStrengthIdx,planeIdx) = trueBeamParams_ref(planeIdx+1,1,1) + planeOffsets(planeIdx);
        end
        
        %Plot side view: True positions
        figure(10)
        for beamIdx = 1:length(offsets)
            plot([0, abs(planePos)], squeeze( trueBeamParams(:,beamIdx,1) )', ...
                'Color', colorCycle(kickStrengthIdx,:))
        end
        h = plot([0, abs(planePos)], squeeze( trueBeamParams_ref(:,1) )', ...
            'Color', colorCycle(kickStrengthIdx,:), 'LineWidth', 4, ...
            'DisplayName', num2str(kickStrengths(kickStrengthIdx)));
        refPlotHandles = [refPlotHandles, h]; %#ok<AGROW>
    end
    
    %% Reconstruct
    
    %Coordinate system definition:
    % - Average position on plane 1 = 0
    % - Average angle between plane 1 and plane 2 = 0
    
    % Average position in plane 1 before first correction;
    % the coordinate system will be shifted so that this becomes 0;
    avgPos0 = mean(hitPositions_ref(:,1))  %#ok<NOPRT> %[mm]
    
    % Average angle in plane 1 (measured between plane 1 and 2) before 2nd correction;
    % the coordinate system will be skewed(~=rotated so that this becomes 0 (and it has been shifted);
    % Note: mean(hitPositions_ref(:,1))-avgPos0 ~= 0; this was included for clarity
    % Note: In a 'reconstructed' case, planePos(1) is most likely 0; here
    avgAng1 = ( (mean(hitPositions_ref(:,2))-avgPos0) - (mean(hitPositions_ref(:,1))-avgPos0)) / (planePos(2)-planePos(1)) %#ok<NOPRT> %[mm/m = mrad]
    
    %Compute "recentered" positions and angles:
    beam_ref = zeros(2,length(kickStrengths));
    beam_ref(1,:) = hitPositions_ref(:,1)-avgPos0; %[mm]
    beam_ref(2,:) = ( (hitPositions_ref(:,2)-avgPos0) - (hitPositions_ref(:,1)-avgPos0) ) / (planePos(2)-planePos(1)) - avgAng1; %[mm/m = mrad]
    beam_ref %#ok<NOPRT>
    
    assert( abs(mean(beam_ref(1,:))) < 1e-14 );
    assert( abs(mean(beam_ref(2,:))) < 1e-14 )
    
    %Where would the beam with offset=0 end up in the structure?
    ballisticBiasStructure = beam_ref(1,:) + beam_ref(2,:)*(abs(planePos(planePos_structureIdx))-planePos(1));
    %Where would the beam with offset=0 end up on the final screen?
    ballisticBiasScreen = beam_ref(1,:) + beam_ref(2,:)*(planePos(end)-planePos(1));
        
    %% Plot
    
    figure(1)
    clf();
    hold on;
    for kickStrengthIdx=1:length(kickStrengths)
        plot( offsets, squeeze(hitPositions(kickStrengthIdx,end,:)), ...
             'DisplayName', num2str(kickStrengths(kickStrengthIdx)) );
    end
    grid on;
    title('Raw data')
    xlabel('Initial offset [mm]')
    ylabel('Position on screen [mm]')
    lgd = legend('Location', 'northwest');
    title(lgd, 'Kick [mrad/mm]')
    print('offset_raw.png', '-dpng')
    
    figure(2)
    clf()
    hold on;
    for kickStrengthIdx=1:length(kickStrengths)
        plot( offsets, squeeze(hitPositions(kickStrengthIdx,end,:))-ballisticBiasScreen(kickStrengthIdx), ...
             'DisplayName', num2str(kickStrengths(kickStrengthIdx)) );
%         plot( offsets, squeeze(hitPositions2(kickStrengthIdx,end,:)), ...
%             'DisplayName', num2str(kickStrengths(kickStrengthIdx)) );
    end
    grid on
    title('Correction 1')
    xlabel('Initial offset [mm]')
    ylabel('Consistent position on screen [mm]')
    lgd = legend('Location', 'northwest');
    title(lgd, 'Kick [mrad/mm]')
    print('offset_corr1.png', '-dpng')

    figure(3)
    clf()
    hold on;
    for kickStrengthIdx=1:length(kickStrengths)
        plot( offsets + ballisticBiasStructure(kickStrengthIdx), ...
              squeeze(hitPositions(kickStrengthIdx,end,:))-ballisticBiasScreen(kickStrengthIdx)+ballisticBiasStructure(kickStrengthIdx), ...
              'DisplayName', num2str(kickStrengths(kickStrengthIdx)) );
    end
    grid on;
    title('Correction 2')
    xlabel('Consistent position in structure [mm]')
    ylabel('Consistent position on screen [mm]')
    lgd = legend('Location', 'northwest');
    title(lgd, 'Kick [mrad/mm]')
    print('offset_corr2.png', '-dpng')
    
    figure(10)
    xlabel('s [m]')
    ylabel('True x [mm]')
    for planeIdx = 1:length(planePos)
        p = planePos(planeIdx);
        if p >= 0
            xline(p, 'k--', 'LineWidth',2);
        else
            assert (planeIdx == planePos_structureIdx)
            xline(-p,'r--', 'LineWidth',2);
            plot([-p-0.2, -p+0.2], [planeOffsets(planeIdx), planeOffsets(planeIdx)], 'r--', 'LineWidth', 2);
        end
        
    end
    title('Side view')
    lgd = legend(refPlotHandles, 'Location', 'northwest');
    title(lgd, 'Kick [mrad/mm]')
    print('sideview.png', '-dpng')

end

