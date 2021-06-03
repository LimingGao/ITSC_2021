function fcn_plotDrivingParameters(error_lookahead,delta,tout)

%
% figure(512347)
% plot(tout,error_hdg*180/pi)
% title('Heading Error')
% xlabel('Time (s)')
% ylabel('Error (deg)')
% ylim([0 2*pi])

figure(512348)
plot(tout,error_lookahead);
title('Lateral Offset Error');
xlabel('Time (s)');
ylabel('Error (m)');

figure(512349)
plot(tout,delta*180/pi);
title('Vehicle Steering Response');
xlabel('Time (s)');
ylabel('Steering Angle (deg)');

end