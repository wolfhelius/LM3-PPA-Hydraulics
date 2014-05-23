clear all
[timestamp_day, hh, ~, ~, psi_rs(:,1), psi_rs(:,2), psi_rs(:,3)] = ...
    textread('../work/WOLF/vegn_cc_psi_rs', '%s %f %f %f %f %f %f', 'headerlines', 2, 'delimiter', ' ,:');

[timestamp_day, hh, ~, ~, psi_r(:,1), psi_r(:,2), psi_r(:,3)] = ...
    textread('../work/WOLF/vegn_cc_psi_r', '%s %f %f %f %f %f %f', 'headerlines', 2, 'delimiter', ' ,:');

[timestamp_day, hh, ~, ~, psi_x(:,1), psi_x(:,2), psi_x(:,3)] = ...
    textread('../work/WOLF/vegn_cc_psi_x', '%s %f %f %f %f %f %f', 'headerlines', 2, 'delimiter', ' ,:');

[timestamp_day, hh, ~, ~, psi_l(:,1), psi_l(:,2), psi_l(:,3)] = ...
    textread('../work/WOLF/vegn_cc_psi_l', '%s %f %f %f %f %f %f', 'headerlines', 2, 'delimiter', ' ,:');

[yyyy, timestamp_day] = strtok(timestamp_day, '-'); yyyy = str2double(yyyy);
[mo, dd] = strtok(timestamp_day, '-'); mo = str2double(mo)-1; dd = -str2double(dd)-1;
doy = yyyy+mo/12+dd/365;
tod = hh/24;
%clear timestamp_day yyyy mo dd hh 
dectime = doy+tod/365;
% figure
% subplot(411);   plot(dectime, psi_rs); title('psi rs, Mpa');
% subplot(412);   plot(dectime, psi_rs); xlim([1950 1951-eps]);
% subplot(413);   plot(dectime, psi_rs); xlim([1950.5 1950.6]);
% subplot(414);   plot(dectime, psi_rs); xlim([1950.53 1950.54]);

figure
subplot(411);   plot(dectime, [psi_rs(:,1),psi_r(:,1),psi_x(:,1),psi_l(:,1)]); title('psi , Mpa');
subplot(412);   plot(dectime, [psi_rs(:,1),psi_r(:,1),psi_x(:,1),psi_l(:,1)]); xlim([1950 1951-eps]);
                legend('Root-soil','After root','After woody xylem','After leaf','Location','West');
subplot(413);   plot(dectime, [psi_rs(:,1),psi_r(:,1),psi_x(:,1),psi_l(:,1)]); xlim([1950.5 1950.6]);
subplot(414);   plot(dectime, [psi_rs(:,1),psi_r(:,1),psi_x(:,1),psi_l(:,1)]); xlim([1950.53 1950.54]);