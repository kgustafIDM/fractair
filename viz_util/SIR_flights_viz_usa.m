% using total_indiv doesn't make sense anymore since individuals can transfer
% disp('total number recovered or still infected up until time t:')
% final_size = repmat(total_indiv,[1,1,num_cities])-final_state(species1,end)
%
cityset = [tophubsUSA+1; infecthub];
%figure;
counter = 1;
for idcity=1:size(cityset,1)
    figure;%subplot(size(cityset,1),1,counter);
    counter = counter + 1;
    plot(time_vector, final_state(species1,1:end,cityset(idcity)), time_vector, final_state(species2,1:end,cityset(idcity)));
    xlabel('Time'); ylabel('count'); % legend('Susceptible','Infected');
end

load /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/withinUSA.mat DepartureAirport DLat DLong
%%
tenpercent = 0.1*total_indiv;

% first find the time of the last occurrence of the global maximum value of infected
for ii = 1:size(final_state,3)
    maxinfect(ii) = max(final_state(2,:,ii));
    maxtimes = find(final_state(2,:,ii)==maxinfect(ii));
    maxhot(ii) = maxtimes(end);
    maxpercent(ii) = max(final_state(2,:,ii))./sum(final_state(:,1,ii),1);
end

% then locate the time at which less than 10 percent is infected, after the
% maximum, so this is on the recovery phase
for ii = 1:size(final_state,3)
    if final_state(2,end,ii)<tenpercent
        downs = find(final_state(2,maxhot(ii):end,ii)<tenpercent);
    else
        downs = size(final_state,2);
    end
    belowtenpercent(ii) = maxhot(ii)+downs(1)-1;
end

% to determine when more than 10 percent of the population is infected
for ii = 1:size(final_state,3)
    ups = find(final_state(2,:,ii)>tenpercent);
    neverabove = isempty(ups);
    if neverabove==1
        abovetenpercent(ii)=size(final_state,2);
    else
        abovetenpercent(ii) = ups(1);
    end
end

% histograms of the critical junctures times for all airport cities
nbins = 20;
[abins,abovehist]=hist(time_vector(abovetenpercent),nbins);
[mbins,infmaxhist]=hist(time_vector(maxhot),nbins);
[bbins,belowhist]=hist(time_vector(belowtenpercent),nbins);

% plot the before, max and after epidemic on the same axes
figure; plot(abovehist,abins,'r'); hold on; plot(infmaxhist,mbins,'b'); plot(belowhist,bbins,'k');
%%
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
range_infected = range(final_state(2,:,:),2);
notmin_rangeinf = find(range_infected>min(range_infected));
clear str1
filename = 'TESTusagfluSIS5_hub44_new.gif';
timeavg = squeeze(mean(final_state(2,:,:),2));
for n = firstframe:lastframe
    if n == firstframe;
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
        scatter(xallIATA(infecthub),yallIATA(infecthub),280,'g+'); 
        scatter(xallIATA(infecthub),yallIATA(infecthub),280,'go');
        str1(1) = {['year ',num2str(floor(mod(time_vector(n),52)))]};
        htlabel = text(610,600,str1,'Color','w','FontSize',30);
        hcbar=colorbar('East'); set(hcbar,'YTickLabel',[]);
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
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
                    scatter(xallIATA(citycheck),yallIATA(citycheck),topdot,colorfactor.*squeeze(final_state(2,n,citycheck)),'filled');
                else
                    scatter(xallIATA(citycheck),yallIATA(citycheck),otherdot,colorfactor.*squeeze(final_state(2,n,citycheck)),'filled');
                end
            end
            % add special emphasis to the initially infected hub
            scatter(xallIATA(infecthub),yallIATA(infecthub),280,'g+'); 
            scatter(xallIATA(infecthub),yallIATA(infecthub),280,'go');
            str1(1) = {['year ',num2str(floor(mod(time_vector(n),52)))]};
%            str1(2) = {['week ',num2str(round(time_vector(n)*52))]};            
            htlabel = text(610,600,str1,'Color','w','FontSize',30);
            hcbar = colorbar('East'); set(hcbar,'YTickLabel',[]);
            drawnow
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'gif','WriteMode','append');
            delete(htlabel);
        end
    end
end
%%

close all
figure;
filename = 'hist_USAhub22.gif';
for n = 1:50:size(final_state,2)
    hist(squeeze(final_state(2,n,:)),[0:1:100]);
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if n == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

figure;
plot(time_vector_SIR, final_state_SIR(species1,1:end,(idcity)),'r',time_vector_SIR, final_state_SIR(species2,1:end,(idcity)),'k');

figure;
plot(time_vector_SIS, final_state_SIS(species1,1:end,(idcity)),'r',time_vector_SIS, final_state_SIS(species2,1:end,(idcity)),'k');

% plot the paths with some interesting colors
% for lk=1:size(xMP,1)
%     scatter(xMP(lk,:),yMP(lk,:),500.*[1:1:size(LonPath,2)],25.*[1:1:size(LonPath,2)],'LineWidth',3)
%     plot(xMP(lk,1:10),yMP(lk,1:10),'r--','LineWidth',2);
%     labels = flightpath(lk,:);
%     text(xMP(lk,:), yMP(lk,:), labels, 'Color','w', 'VerticalAlign','bottom', 'HorizontalAlign','right','FontSize',18)
% end