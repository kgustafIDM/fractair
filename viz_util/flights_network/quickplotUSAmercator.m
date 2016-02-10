
% find latitude and longitude for all airports
[uniqAirport,idxuniq] = unique(DepartureAirport);
uniqLat = DLat(idxuniq);
uniqLong = DLong(idxuniq);

% Code from http://stackoverflow.com/users/97160/amro
% world map in Mercator projection
fnameMerc = '/Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/Mercator-projectionlowsat.jpg';
imgMerc = imread(fnameMerc);
[imgHMerc,imgWMerc,~] = size(imgMerc);

% Mercator projection of a path
%[xMP,yMP] = mercatorProjection(LonPath, LatPath, imgWMerc, imgHMerc);

% Mercator projection of all airport in USA
[xallIATA,yallIATA] = mercatorProjection(uniqLong,uniqLat,imgWMerc,imgHMerc);

colormap jet
% plot the entire globe projected on a Mercator map
imshow(imgMerc, 'InitialMag',45,'Border','tight');
hold on;
% zoom in on InDia
axis([1370 1600 570 760])
% zoom in on USA
axis([300 675 450 675])

% size of top hub dots and size of all other dots
topdot = 600; otherdot = 300;
% factor for scaling up the colorbar
colorfactor = 350;

% movie timing parameters
firstframe = 1;
framestep = 10;
lastframe = size(final_state,2);
studtime = time_vector(firstframe);

% secgrowth_inf = zeros(size(final_state,3),2);
% for j = 1:size(final_state,3)
%     secgrowth_inf(j,:) = polyfit(time_vector,final_state(2,:,j),1);
% end
% find the airports that show no transfers of passengers with the reduced
% proxy populations
range_infected = range(final_state(2,:,:),2);
notmin_rangeinf = find(range_infected>min(range_infected));
%
clear str1
timeavg = squeeze(mean(final_state(2,:,:),2));
for n = firstframe:lastframe
    
    if (time_vector(n) > studtime + 6/52)
        n
        time_vector(n)
        studtime = time_vector(n);
        for city=1:size(notmin_rangeinf,1)
            citycheck = notmin_rangeinf(city);
            % make the color proportional to the infection
            if intersect(citycheck,tophubsUSA)
                scatter(xallIATA(citycheck),yallIATA(citycheck),topdot,colorfactor.*squeeze(final_state(2,n,citycheck)),'filled');
            else
                scatter(xallIATA(citycheck),yallIATA(citycheck),otherdot,colorfactor.*squeeze(final_state(2,n,citycheck)),'filled');
            end
        end
        % add special emphasis to the initially infected hub
        year = floor(mod(time_vector(n),52))+1;
        scatter(xallIATA(tophubsUSA(year)),yallIATA(tophubsUSA(year)),280,'g+');
        scatter(xallIATA(tophubsUSA(year)),yallIATA(tophubsUSA(year)),280,'go');
        str1(1) = {['year ',num2str(year)]};
        %            str1(2) = {['week ',num2str(round(time_vector(n)*52))]};
        htlabel = text(610,600,str1,'Color','w','FontSize',30);
        hcbar = colorbar('East'); set(hcbar,'YTickLabel',[]);
        drawnow
        delete(htlabel);
    end
end