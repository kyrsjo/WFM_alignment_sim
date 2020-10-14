function [trueBeamParams] = kickSim(initialBeamParams,kickStrength, planePos, verbose, structureTransOffset)
    % kickSim: Simulate transport of beams through screens and kicks
    %
    % Input:
    % - initialBeamParams : List of [ [x,xp]' , ...] initial beam parameters [mm, mrad]
    % - kickStrength      : List of kick strengths in [mrad/mm] for each beam
    % - planePos: Array of plane positions [m].
    %             Negative => accel structure that gives kick (at position abs(planePos) )
    %             Positive : observation point
    % - verbose: If present and true, output debug info
    % - structureTransOffset: Transverse offset of structure center [mm]
    %            relative to "true" coordinate system
    % Output:
    % - trueBeamParams:    Positions and angles of each beam in each plane,
    %                      incluing the initial position and angle [mm,mrad]
    %                      (size 2 x (number of planes+1), first index = initial position )
    
    assert ( size(initialBeamParams,2) == length(kickStrength) )
    
    if ~exist('verbose','var')
       verbose = false; 
    end
    if ~exist('structureTransOffset','var')
        structureTransOffset = 0;
    end
    
    if verbose
        disp('initialBeamParams [mm; mrad]:')
        disp(initialBeamParams)
        disp('kickStrength: [mrad/mm]')
        disp(kickStrength)
        disp('planePos [m]:')
        disp(planePos)    
    end
    
    %True beam parameters in each plane;
    % indices: planeIdx , beamIdx, x/xp
    trueBeamParams = zeros(length(planePos)+1,size(initialBeamParams,2),2);
    % Plane 1 is always the initial parameters
    planePos = [0.0, planePos];
    for beamIdx=1:size(initialBeamParams,2)
        trueBeamParams(1,beamIdx,:) = initialBeamParams(:,beamIdx);
    end
    
    %Compute the true beam parameters
    for planeIdx = 2:length(planePos)
        L = abs(planePos(planeIdx)) - abs(planePos(planeIdx-1)); %Drift length to this plane[m]

        for beamIdx = 1:size(initialBeamParams,2)
            %drift to current plane
            trueBeamParams(planeIdx, beamIdx, 1)     = trueBeamParams(planeIdx-1, beamIdx, 1) + trueBeamParams(planeIdx-1, beamIdx, 2)*L;
            trueBeamParams(planeIdx, beamIdx, 2)     = trueBeamParams(planeIdx-1, beamIdx, 2);
            %kick!
            if planePos(planeIdx) < 0
                trueBeamParams(planeIdx, beamIdx, 2) = trueBeamParams(planeIdx,beamIdx, 2) + ...
                    (trueBeamParams(planeIdx,beamIdx, 1) - structureTransOffset)*kickStrength(beamIdx);
            end
        end
    end

    if verbose
        disp('Final beam params (planeIDx,beamIDx,x/xp):')
        disp(trueBeamParams)
    end

end