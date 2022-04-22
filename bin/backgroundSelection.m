function backgroundSelection(Ntraces, len, donor, acceptor)

%{



%}





    % define FRET
    fret=zeros(Ntraces/2,len);
    for j=1:Ntraces/2
        for k=1:len
            if donor(j,k)+acceptor(j,k)==0
                fret(j,k)=0.5;
            else fret(j,k)=gamma.*(acceptor(j,k)-leakage*donor(j,k))/(donor(j,k)+gamma.*acceptor(j,k));
            end
            if fret(j,k)>1.5
                fret(j,k)=1.5;
            end
            if fret(j,k)<-0.5
                fret(j,k)=-0.5;
            end
        end
    end
    
    hdl=gcf;        
    set(hdl,'Position',[4,300,1276,648]);
    
%%%%% JUMP TO SPECIFIC MOLECULE
%%%%% j=250;

j=0;
    while(j<Ntraces/2)
        j=j+1;
        
        % plot donor, acceptor, and FRET
        h1=subplot(5,1,[1,2,3]);
        plot(time,donor(j,:),'g',time,(acceptor(j,:))-abs(leakage*donor(j,:)),'r', time,(acceptor(j,:)+donor(j,:)+200), 'k' );
        ylim([-25 max((acceptor(j,:)+donor(j,:)+400))]);
        xlim([-20*timeunit max(time)+100*timeunit]);
        
        title(['File:' fname '   Molecule:' int2str(j-1)]);
        grid on;
        zoom on;
        h2=subplot(2,1,2);
        plot(time,fret(j,:),'b');
         ylim([-0.02 1.02]);
         xlim([-20*timeunit max(time)+100*timeunit]);
         
        grid on;
        zoom on;
        linkaxes([h1,h2], 'x');
        %frames_to_avg=3;
        
        option=input('skip (enter), Work on this trace(s), go back (b), save trace (t) ','s');
        if option=='t'
            fname1=[fname ' tr' num2str(j-1) '.dat'];
            output=[time' donor(j,:)' acceptor(j,:)'];
            save(fname1,'output','-ascii') ;
            j=j-1;
        else
        
%         if option=='s'
%             disp('Molecule skipped');
%         else
            if option=='b'
                if j==1
                    disp('This is the first molecule');
                    j=j-1;
                else
                    disp('Going to previous molecule');
                    j=j-2;
                end
            else
                if option=='c'
                    if (j==1 && i==3)
                        disp('This is the first molecule');
                    else
                        if number_of_frames==0
                            disp('There is no selection to delete');
                        else
                            if was_deleted==0
                                total_N_of_frames=total_N_of_frames-number_of_frames;                        
                                disp('Previous selection was cancelled');
                                was_deleted=1;
                            else
                                disp('Previous selection was already cancelled');
                            end
                        end
                    end
                    j=j-1;
                else
                    if option=='s'
                    disp('Click twice to select the range (press enter to skip)');
                    [x,y]=ginput(2);
                    
                    if isempty(x) | size(x)==1
                        x=ones(2);
                    else
                        x(1)=round(x(1)/timeunit);
                        x(2)=round(x(2)/timeunit);
                        if x(1)<1
                            x(1)=1;
                        end
                        if x(2)<1
                            x(2)=1;
                        end
                        if x(1)>len
                            x(1)=len;
                        end
                        if x(2)>len
                            x(2)=len;
                        end
                    end
                    
                    if x(1)>x(2)
                        temp=x(1);
                        x(1)=x(2);
                        x(2)=temp;
                    end
                    
                    disp(['Selected range is frame ',num2str(x(1)),' to frame ',num2str(x(2))]);
                    number_of_frames=x(2)-x(1)+1;
                    was_deleted=0;
                    
                     %%%%%%%%%%%%%%
                    x(2)=x(2)-mod((x(2)-x(1)),sm);
                    %%%%%%%%%%%%%
                    

                    % Manual background selection
                    
                    if bg_definition
                        disp('Click 4X  to select the range to calculate background (enter to skip)');
                        [xbc3,ybc3]=ginput(2);
                        [xbc5,ybc5]=ginput(2);
                       
                        if xbc3(2)> xbc3(1);
                            xbc3(1)=round(xbc3(1)/timeunit);
                            xbc3(2)=round(xbc3(2)/timeunit);
                            bg_donor=mean(donor(j,xbc3(1):xbc3(2)));
                        else
                            bg_donor=input('donor background? ');
                        end
                        
                         if xbc5(2)> xbc5(1);
                            xbc5(1)=round(xbc5(1)/timeunit);
                            xbc5(2)=round(xbc5(2)/timeunit);
                            bg_acceptor=mean(acceptor(j,xbc5(1):xbc5(2)) - leakage.*donor(j,xbc5(1):xbc5(2)));
                        else
                            bg_acceptor=input('acceptor background? ');
                        end
                    
                        
                        
                      
                    end
                    
                    disp(['Backgrounds of donor and acceptor are ',num2str(bg_donor),' and ',num2str(bg_acceptor)]);
                    
                     %%%% GAMMA FACTOR CORRECTION  %%%%%%%%%%%%%
                   %%%%%% Either manually enter GAMMA or click 4 times (donor low, donor high then acceptor low, acceptor high)
                   %%%%%%%%%%% let the program measure Gammma. If you don't know Gamma put 0.%%%
                   
                   gamma=0;deltaA=0;deltaD=0; 
                   gc=input('Manual (enter) or automatic GAMMA factor? ', 's');
                    if isempty(gc)
                         gamma=input('What is the gamma factor? ');
                         if isempty (gamma)
                             gamma=0;
                         end;
                    else
                        disp('Click 4X to select the DONOR range to calculate gamma '); 
                        [xdon1,ydon1]=ginput(2);
                        [xdon2,ydon2]=ginput(2);
                        xdon1(1)=round(xdon1(1)/timeunit);  xdon1(2)=round(xdon1(2)/timeunit);
                        xdon2(1)=round(xdon2(1)/timeunit);  xdon2(2)=round(xdon2(2)/timeunit);
                        D1=mean(donor(j,min(xdon1(1),xdon1(2)): max(xdon1(1),xdon1(2))));
                        D2=mean(donor(j,min(xdon2(1),xdon2(2)): max(xdon2(1),xdon2(2))));
                        deltaD=abs(D1-D2);
                        
                        
               
                        disp('Click 4X to select the ACCEPTOR range to calculate gamma '); 
                        [xacc1,yacc1]=ginput(2);
                        [xacc2,yacc2]=ginput(2);
                        xacc1(1)=round(xacc1(1)/timeunit);  xacc1(2)=round(xacc1(2)/timeunit);
                        xacc2(1)=round(xacc2(1)/timeunit);  xacc2(2)=round(xacc2(2)/timeunit);
                        acceptorleak=acceptor-abs(leakage*donor); %leakage corrected acceptor
                        A1=mean(acceptorleak(j,min(xacc1(1),xacc1(2)): max(xacc1(1),xacc1(2))));
                        A2=mean(acceptorleak(j,min(xacc2(1),xacc2(2)): max(xacc2(1),xacc2(2))));
                        deltaA=abs(A1-A2);
                          
                        gamma=deltaA/deltaD;
                        disp(['Gamma factor is ',num2str(gamma)]);
                    end
                     fid = fopen(['gamma' fname '.dat'],'a+');
                     fprintf(fid,'%4.3f\n', gamma);
                     fclose(fid);
                    %%%% END of GAMMA FACTOR CORRECTION%%%%
                    %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if number_of_frames>1
                        
%                         for k=0:(x(2)-x(1)-sm)/sm;
%                             temp_donorsm=mean(donor(j,(x(1)+sm*k):(x(1)+sm*(k+1))))-bg_donor;
%                             temp_acceptorsm=mean(acceptor(j,(x(1)+sm*k):(x(1)+sm*(k+1))))-bg_acceptor;
%                             allfretsm(k+1)=gamma.*(temp_acceptorsm-leakage*temp_donorsm)/(gamma.*temp_acceptorsm+(temp_donorsm));
% %                             saved_frames(1,total_N_of_frames+k)=(temp_acceptor-leakage*temp_donor)/(temp_acceptor+temp_donor);
%                         end
                        temp_donorsm=smooth(donor(j,x(1):x(2)),3)-bg_donor;
                        temp_acceptorsm=smooth(acceptor(j,x(1):x(2)),3)-bg_acceptor;
                        allfretsm=(temp_acceptorsm-abs(leakage.*temp_donorsm))./(temp_acceptorsm+temp_donorsm);
                        [nn,xout] =hist(allfretsm,gap);
                        nn_norm=nn/sum(nn);
                        N_His=(N_His+nn_norm)/sum(N_His+nn_norm);
                        N_His_raw=N_His_raw+nn;
                        %%%%%%%%%%%%%%
                         output=[xout' , N_His', N_His_raw'];
                         time_fname=['selectedFRET_SMave' fname '.dat'];
                         save(time_fname,'output','-ascii') ;

                        
                         %%%%%%%%%%%%%%%%%%SM average.
%                         allfret=saved_frames(1,total_N_of_frames+1:total_N_of_frames+number_of_frames);
%                         fid = fopen(['selectedFRET_SMave' fname '.dat'],'a+');
%                         fprintf(fid, '%4.4f\n', allfretsm);
%                          fclose(fid);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        mean_donor=mean(donor(j,x(1):x(2)))-bg_donor;
                        mean_acceptor=mean(acceptor(j,x(1):x(2)))-bg_acceptor;
                        allfret=(mean_acceptor-abs(leakage*mean_donor))/(mean_acceptor+mean_donor)
                        allfret_gamma=(mean_acceptor-abs(leakage*mean_donor))/(mean_acceptor+gamma*mean_donor)
                        selfret=(acceptor(j,x(1):x(2))-bg_acceptor-abs(leakage*(donor(j,x(1):x(2)))))./((donor(j,x(1):x(2)))+acceptor(j,x(1):x(2))-bg_donor-bg_acceptor);
                        seltotal=donor(j,x(1):x(2))+acceptor(j,x(1):x(2));
                        totalI=(mean_acceptor+mean_donor)
                        %%%%%%%%%%%%%%%%%%I just added these 4 lines.
%                         allfret=saved_frames(1,total_N_of_frames+1:total_N_of_frames+number_of_frames);
                        fid = fopen(['selectedFRETavseg' fname '.dat'],'a+');
                        fprintf(fid, '%4.0f\t %5.1f\t  %4.3f\t %4.3f\n', j,totalI,allfret,gamma);
                         fclose(fid);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        selres=[seltotal ; selfret];
                        fid = fopen(['selectedFRET' fname '.dat'],'a+');
                        fprintf(fid,'%5.1f  %4.3f\n', selres);
                        fclose(fid);
                        
                         
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %%%This is the selected traces that we will use for
                       %%%future analysis. the columns are: TIME, DONOR
                       %%%(background corrected), Acceptor (background and
                       %%%leakage corrected) and gamma*DONOR. In case of
                       %%%Gamma, FRET=A/(A+gamma*D).
                    
                        fname5=[fname ' tr' num2str(j-1) '.dat'];
                        output5=[time(x(1):x(2))'  (donor(j,x(1):x(2))-bg_donor)' (acceptor(j,x(1):x(2))-abs(leakage*(donor(j,x(1):x(2))))-bg_acceptor)' gamma*(donor(j,x(1):x(2))-bg_donor)'];
                        save(fname5,'output5','-ascii') ;
                        
%                          fname1=[fname ' tr' num2str(j) '.dat'];
%                          output=[time(x(1):x(2))' (donor(j,x(1):x(2))-bg_donor)' gamma.*((acceptor(j,x(1):x(2))-bg_acceptor)-leakage*(donor(j,x(1):x(2))-bg_donor))'];
%                          save(fname1,'output','-ascii') ;
%                          plot(time,(donor(j,:)-bg_donor),'g',time,gamma.*((acceptor(j,:)-bg_acceptor)-leakage*(donor(j,:)-bg_donor)),'r', time,(gamma.*acceptor(j,:)+donor(j,:)+1000), 'k' );

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         figure;hist(allfret)
                        total_N_of_frames=total_N_of_frames+number_of_frames;
                        allfret=[];
                        
                        more_range=input('Press b to select more range [default=no] ','s');
                        if more_range=='b'
                            j=j-1;
                        end
                    end
                    
                end
            end
            end
        end
    end