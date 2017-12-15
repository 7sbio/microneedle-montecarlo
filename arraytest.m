%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo Simulation of Microneedle Insertions
%
% (c) 2017 Seventh Sense Biosystems Inc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define some constants
N = 100;                              % Number of configurations to try
needleN = 30;                         % Number of needles in the array
needleXspacing = 0.2;                 % Spacing between needles in the x-direction in mm
needleYspacing = 0.85;                % Spacing between needles in the y-direction in mm
needleXlength = 0.350;                % Size of the needles in the x-dimension in mm
needleYlength = 0.05;                 % Size of the needles in the y-dimension in mm
needleRows = 6;                       % Number of rows of needles
capillaryDensity = 14;                % Number of capillaries per square mm of skin
capillaryDiameter = 0.015;            % Diameter of a capillary loop in mm
maxiter = 500;                        % Maximum number of iterations before throwing out a trial
capillaryMinimumDistance = 0.080;     % Minimum inter-capillary distance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some parameters based on defined constants
boxXlength = (needleN/needleRows)*(needleXspacing+needleXlength);
boxYlength = needleRows*(needleYspacing+needleYlength);
if boxXlength < 1, boxXlength = 1.; end
if boxYlength < 1, boxYlength = 1.; end
  
%   
needleCenterXvalues = zeros(1,needleN);
needleCenterYvalues = zeros(1,needleN);
for i=1:needleN,
    needleCenterXvalues(i) = (i-(ceil(i/(needleN/needleRows))-1)*(needleN/needleRows))*(needleXspacing+needleXlength)-0.5*(needleXspacing+needleXlength);
    needleCenterYvalues(i) = ceil(i/(needleN/needleRows))*(needleYspacing+needleYlength)-0.5*(needleYspacing+needleYlength);
end
capillaryN = ceil(boxXlength*boxYlength*capillaryDensity);

% Arrays for histogram of capillary strikes
direct = zeros(1,capillaryN);
glance = zeros(1,capillaryN);
total = zeros(1,capillaryN);


halfBoxXlength = 0.5*boxXlength;
halfBoxYlength = 0.5*boxYlength;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The main loop of the program
tic
for i = 1:N;
    % Display trial number 
    disp(i);fflush(stdout);
    
    % Initialize arrays for locations of capillaries
    capillaryCenterXvalues = zeros(1,capillaryN);
    capillaryCenterYvalues = zeros(1,capillaryN);
    
    % Place the capillaries
    bad = 0; % flag to set if number of iterations reaches maxiter
    for j = 1:capillaryN
        done = 0; % flag to set when new capillary location is acceptable
        ctr = 0; % counter to keep track of number of iterations used
        while done == 0 && ctr < maxiter
            ctr = ctr+1;
            new_capillaryCenterXvalue = rand*boxXlength; % place capillary randomly in box
            new_capillaryCenterYvalue = rand*boxYlength;
            done = 1; % accept capillary location unless too close to another 
            
            % Check all pairwise capillary-capillary distances 
            if j > 1
                for k = 1:j-1
                    xdist = abs(new_capillaryCenterXvalue-capillaryCenterXvalues(k));
                    ydist = abs(new_capillaryCenterYvalue-capillaryCenterYvalues(k));
                    if xdist > halfBoxXlength, xdist = xdist-boxXlength; end % wrap capillary pattern around box
                    if ydist > halfBoxYlength, ydist = ydist-boxYlength; end
                    dist = sqrt(xdist*xdist+ydist*ydist);
                    if dist < capillaryMinimumDistance
                      done = 0; % if one distance too short, don't accept location
                      continue; % if one distance too short, don't check the rest  
                    end
                end
            end
        end
        capillaryCenterXvalues(j) = new_capillaryCenterXvalue;
        capillaryCenterYvalues(j) = new_capillaryCenterYvalue;
        % If we're at maxiter iterations and still not done, give up
        if ctr == maxiter, 
            bad = 1;
            continue;
        end
    end
    
    if bad==1, continue; end; % If we couldn't generate capillary pattern, don't count intersections

    % Count the intersections between needles and capillaries
    trialDirectStrikes = 0;
    trialGlanceStrikes = 0;
    trialTotalStrikes = 0;
    for j = 1:capillaryN,
        for k = 1:needleN;
            xdist = abs(needleCenterXvalues(k)-capillaryCenterXvalues(j));
            ydist = abs(needleCenterYvalues(k)-capillaryCenterYvalues(j));
            
            if xdist > halfBoxXlength, xdist = abs(xdist-boxXlength); end % wrap capillary pattern around box
            if ydist > halfBoxYlength, ydist = abs(ydist-boxYlength); end % wrap capillary pattern around box
            
            if xdist < 0.5*(needleXlength+capillaryDiameter) && ydist < 0.5*(needleYlength+capillaryDiameter),
                trialTotalStrikes = trialTotalStrikes+1;
                if xdist < 0.5*(needleXlength-capillaryDiameter) && ydist < 0.5*(needleYlength-capillaryDiameter),
                  trialDirectStrikes = trialDirectStrikes+1;
                else
                  trialGlanceStrikes = trialGlanceStrikes+1;
                end
            end
        end                    
    end

    % If trial is good, add strike counts to histogram bins
    direct(trialDirectStrikes+1) = direct(trialDirectStrikes+1)+1;
    glance(trialGlanceStrikes+1) = glance(trialGlanceStrikes+1)+1;
    total(trialTotalStrikes+1) = total(trialTotalStrikes+1)+1;
end

% Write the output to a file
str = strcat(int2str(needleN),'MN_',int2str(capillaryDensity),'cap_per_mm2.csv');
fid1 = fopen(str,'w');
fprintf(fid1,'N = %g\n',N);
fprintf(fid1,'needleN = %g\n',needleN);
fprintf(fid1,'needleXspacing = %g\n',needleXspacing);
fprintf(fid1,'needleYspacing = %g\n',needleYspacing);
fprintf(fid1,'needleXlength = %g\n',needleXlength);
fprintf(fid1,'needleYlength = %g\n',needleYlength);
fprintf(fid1,'needleRows = %g\n',needleRows);
fprintf(fid1,'capillaryDensity = %g\n',capillaryDensity);
fprintf(fid1,'capillaryDiameter = %g\n',capillaryDiameter);
fprintf(fid1,'maxiter = %g\n',maxiter);
fprintf(fid1,'capillaryMinimumDistance = %g\n',capillaryMinimumDistance);

fprintf(fid1,'Capillaries,Direct,Glance,Total\n');
for i=1:capillaryN,
    fprintf(fid1,'%6.4d,%6.4d,%6.4d,%6.4d\n',i-1,direct(i),glance(i),total(i));
end
fclose(fid1);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
