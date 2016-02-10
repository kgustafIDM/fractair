clear str1
filename = 'usanew43scans.gif';

Drange = [2 5 10 20 50 80 100];
alpharange = [0.1:0.2:0.9 1 1.1:0.2:2];
end_times = 0.05:0.05:10;

load ~/Research/EMOD/SISfrac/SISUSA_4_3_nocycle.mat final_state time_vector

cd ~/Desktop/paperBayati/results/eucdisterror/newparamscans/

% movie timing parameters
firstframe = 1;
%framestep = 10;
lastframe = size(end_times,2);
%

f2 = figure('Visible','on');
axes1 = axes('Parent',f2,'FontSize',20);
%%
for n = firstframe:lastframe
    endtimeid = find(min(abs(time_vector-end_times(n)))==abs(time_vector-end_times(n)));
    sigUSA43 = squeeze(mean(final_state(2,1:endtimeid,:),3))./1000;
    sizesigUSA43 = size(sigUSA43,2);
    array = zeros(size(Drange,2),size(alpharange,2));
    l = 0;
    for Dalpha = Drange
        l = l+1;
        k = 0;
        for alphafrac = alpharange
            k = k+1;
            fname = sprintf('frac43_%.2f_%.2f.mat',Dalpha,alphafrac);
            load(fname,'sigwork');
            sizesigwork = size(sigwork,2);
            if sizesigUSA43<sizesigwork
                resamp = floor(linspace(1,sizesigwork,sizesigUSA43));
                sigwork = sigwork(1,resamp);
            elseif sizesigUSA43>sizesigwork
                resamp = floor(linspace(1,sizesigUSA43,sizesigwork));
                sigUSA43 = sigUSA43(1,resamp);
            else
                sigUSA43 = sigUSA43;
                sigwork = sigwork;
            end
            eucdistUSA = pdist([sigUSA43;sigwork]);
            array(l,k) = eucdistUSA;
        end
    end
    if n == firstframe;
        surf(log10(Drange),alpharange,array');
        axis([0 3 0 2 0 20]);
        str1(1) = {['time ',num2str(end_times(n))]};
        title(str1); xlabel('log_10D'); ylabel('\alpha'); zlabel('L2 norm');
        set(gca,'ZScale','log');
        view(axes1,[34.5 42]); sdf('font20');
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % this initial call of imwrite is necessary for the gif format
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf)
    else
        surf(log10(Drange),alpharange,array');
        str1(1) = {['time ',num2str(end_times(n))]};
        title(str1); xlabel('log_10D'); ylabel('\alpha'); zlabel('L2 norm');
        axis([0 3 0 2 0 20]);
        set(gca,'ZScale','log');
        view(axes1,[34.5 42]); sdf('font20'); 
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end