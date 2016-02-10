load /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/withinUSA.mat

%%
addpath ~/Research/Naef/paperHiSeq/matlabscripts/
load /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/transferUSAnodes.mat

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
[xnodeIATA,ynodeIATA] = mercatorProjection(uniqLong(idnodes(nodeidUSA)),uniqLat(idnodes(nodeidUSA)),imgWMerc,imgHMerc);

colormap jet
% plot the entire globe projected on a Mercator map
imshow(imgMerc, 'InitialMag',45,'Border','tight');
hold on;
% zoom in on InDia
axis([1370 1600 570 760])
% zoom in on USA
axis([300 675 450 675])

% size of top hub dots and size of all other dots
topdot = 500; otherdot = 200;
% factor for scaling up the colorbar
colorfactor = 1;%350;

% if pulling data from smaller, individual node-seeding parameter scan files then let
final_state = finals; final_state = final_state;
time_vector = timev;
% AND change all final_state(2, to final_state(1,

% movie timing parameters
firstframe = 1;
framestep = 1;
lastframe = size(final_state,2);
studtime = time_vector(firstframe);

% secgrowth_inf = zeros(size(final_state,3),2);
% for j = 1:size(final_state,3)
%     secgrowth_inf(j,:) = polyfit(time_vector,final_state(2,:,j),1);
% end
% find the airports that show no transfers of passengers with the reduced
% proxy populations

range_infected = range(final_state(1,:,:),2);
notmin_rangeinf = find(range_infected>min(range_infected));

minmean = min(mean(squeeze(final_state(1,:,:)),2));
maxmean = max(mean(squeeze(final_state(1,:,:)),2));
%
clear str1
filename = 'SISnoisycycleUSA_6_3.gif';
timeavg = squeeze(mean(final_state(1,:,:),2));
for n = firstframe:framestep:lastframe
    if n == firstframe;
        for city=1:size(notmin_rangeinf,1)
            citycheck = notmin_rangeinf(city);
            % make the color proportional to the infection
            if intersect(citycheck,tophubsUSA)
                scatter(xnodeIATA(citycheck),ynodeIATA(citycheck),topdot,colorfactor.*(squeeze(final_state(1,n,citycheck))),'filled');
            else
                scatter(xnodeIATA(citycheck),ynodeIATA(citycheck),otherdot,colorfactor.*(squeeze(final_state(1,n,citycheck))),'LineWidth',3);
            end
        end
        % add special emphasis to the initially infected hub
        scatter(xnodeIATA(infecthub),ynodeIATA(infecthub),280,'k+');
        scatter(xnodeIATA(infecthub),ynodeIATA(infecthub),280,'ko');
        %        str1(1) = {['year ',num2str(floor(mod(time_vector(n),1)))]};
        str1(1) = {['week ',num2str(round(time_vector(n)*52))]};
        htlabel = text(610,600,str1,'Color','w','FontSize',30);
        hcbar=colorbar('East'); ylabel(hcbar,'log2(I_i)'); caxis([0 max(squeeze(final_state(1,1,:)))]); %set(hcbar,'YTickLabel',[]); 
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % this initial call of imwrite is necessary for the gif format
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf)
        delete(htlabel);
    else
        if (time_vector(n) > studtime + 4/52)
            n
            time_vector(n)
            studtime = time_vector(n);
            for city=1:size(notmin_rangeinf,1)
                citycheck = notmin_rangeinf(city);
                % make the color proportional to the infection
                if intersect(citycheck,tophubsUSA)
                    scatter(xnodeIATA(citycheck),ynodeIATA(citycheck),topdot,colorfactor.*(squeeze(final_state(1,n,citycheck))),'filled');
                else
                    scatter(xnodeIATA(citycheck),ynodeIATA(citycheck),otherdot,colorfactor.*(squeeze(final_state(1,n,citycheck))),'LineWidth',3);
                end
            end
            hold on;
            % add special emphasis to the initially infected hub
            scatter(xnodeIATA(infecthub),ynodeIATA(infecthub),280,'k+');
            scatter(xnodeIATA(infecthub),ynodeIATA(infecthub),280,'ko');
            %            str1(1) = {['year ',num2str(floor(mod(time_vector(n),1)))]};
            str1(1) = {['week ',num2str(round(time_vector(n)*52))]};
            htlabel = text(610,600,str1,'Color','w','FontSize',30);
            hcbar=colorbar('East'); caxis([0 max(squeeze(final_state(1,n,:)))]); %set(hcbar,'YTickLabel',[]);
            drawnow
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'gif','WriteMode','append');
            delete(htlabel);
        end
    end
end