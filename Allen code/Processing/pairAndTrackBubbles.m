res_assumption = 10; % assumption of how good the localization is (um)
maxDetectSpeed = (res_assumption / 1e6) / (1 / P.frameRate); % maximum bubble speed we can detect
maxDistPerFrame = (res_assumption / 1e6); % maximum travel distance between frames (m)