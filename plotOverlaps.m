function plotOverlaps(initialBeamParams0,kickStrengths, offsets, planePos, structOff)
    %PLOTOVERLAPS Plot "position scan" plots with variable kickStrengths
    % - initialBeamParams0: List of initial beam parameters for each kick strength
    % - offsets: offset scan array [mm]
    %
    % Assumption: Only a single structure (i.e. negative planePos)!
    %

    if ~exist('structOff', 'var')
        structOff = 0.0;
    end
    
    %Sanity checks
    assert( length(initialBeamParams0) == length(kickStrengths) )
    
    planePos_structureIdx = -1;
    for i=1:length(planePos)
        if planePos(i) < 0
            assert (planePos_structureIdx==-1) %There can only be one!
            planePos_structureIdx = i;
        end
    end
    
    figure(1); clf()
    figure(2); clf()
    figure(3); clf()
    
    for kickStrengthIdx=1:length(kickStrengths)
        initialBeamParams = ones(2,length(offsets)) .* initialBeamParams0(:,kickStrengthIdx);
        initialBeamParams(1,:) = initialBeamParams(1,:)+offsets;
        
        kickStrength = ones(1,length(offsets)) .* kickStrengths(kickStrengthIdx);
        
        %The scan
        trueBeamParams = kickSim(initialBeamParams,kickStrength, planePos);
        %Reference beam
        %trueBeamParams_ref = kickSim(initialBeamParams0(:,kickStrengthIdx),kickStrength(kickStrengthIdx), planePos)
        
        %Where would the beam with offset=0 end up in the structure?
        ballisticBiasStructure = initialBeamParams0(1,kickStrengthIdx) + initialBeamParams0(2,kickStrengthIdx)*abs(planePos(planePos_structureIdx)) + structOff;
        %Where would the beam with offset=0 end up on the final screen?
        ballisticBiasScreen = initialBeamParams0(1,kickStrengthIdx) + initialBeamParams0(2,kickStrengthIdx)*planePos(end);
        
        %plot position on final screen as function of offset
        figure(1)
        plot(offsets, trueBeamParams(end,:,1), 'DisplayName', num2str(kickStrengths(kickStrengthIdx)) );
        hold on
        grid on
        
        %plot position on final screen, corrected for ballistic bias,
        % as a function of offset
        %If no kick, this brings lines back on top of each other
        figure(2)
        plot(offsets, trueBeamParams(end,:,1) - ballisticBiasScreen, ...
            'DisplayName', num2str(kickStrengths(kickStrengthIdx)) );
        
        hold on
        grid on
        
        %plot position on final screen, corrected for ballistic bias,
        % as a function of offset in structure
        figure(3)
        plot(offsets + ballisticBiasStructure, trueBeamParams(end,:,1)-ballisticBiasScreen+ballisticBiasStructure, ...
            'DisplayName', num2str(kickStrengths(kickStrengthIdx)) );
        hold on
        grid on
        
    end
    
    figure(1)
    title('Raw data')
    xlabel('Initial offset [mm]')
    ylabel('Position on screen [mm]')
    lgd = legend('Location', 'northwest');
    title(lgd, 'Kick [mrad/mm]')
    print('raw.png', '-dpng')
    
    figure(2)
    title('Correction 1')
    xlabel('Initial offset [mm]')
    ylabel('Position on screen [mm] - ballistic bias to screen')
    lgd = legend('Location', 'northwest');
    title(lgd, 'Kick [mrad/mm]')
    print('corr1.png', '-dpng')

    figure(3)
    title('Correction 2')
    xlabel('Offset in structure [mm]')
    ylabel('Position on screen relative to structure (??) [mm] - ballistic bias')
    lgd = legend('Location', 'northwest');
    title(lgd, 'Kick [mrad/mm]')
    print('corr2.png', '-dpng')

end

