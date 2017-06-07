if (0)
    figure
    plot(tt_*1e12,mTM_(:,3),tt_*1e12,mRE_(:,3))
    xlabel('t(ps)');ylabel('mz')
    legend('TM','RE')
end
if (0)%m dynamics of individual atoms
    figure
    plot(tt_*1e12,mTM_(:,1),tt_*1e12,mTM_(:,2),tt_*1e12,mTM_(:,3),'linewidth',2)
    legend('mx','my','mz');
    xlabel('t(ps)');ylabel('mTM')
    figure;
    plot(tt_*1e12,mRE_(:,1),tt_*1e12,mRE_(:,2),tt_*1e12,mRE_(:,3),'linewidth',2)
    legend('mx','my','mz');
    xlabel('t(ps)');ylabel('mRE')
    %legend('TM','RE')
end
if (0)%angle b/w two atoms
    figure;
    plot(tt_(1:end-1)*1e12,mangle_(1:end-1)/pi*180,'linewidth',3);
    xlabel('t(ps)');ylabel('angle')
end
if (0)
    mTMampi=sqrt(mTM_(:,1).^2+mTM_(:,2).^2+mTM_(:,3).^2);
    mREampi=sqrt(mRE_(:,1).^2+mRE_(:,2).^2+mRE_(:,3).^2);
    figure;
    plot(tt_*1e12,mREampi,'linewidth',3)
    xlabel('t(ps)');ylabel('norm(mRE)')
    
end
if (0)%draw M-H loop
    machineselc=1;%\0 room 1\lab
    sweepsele=4;%\1alp \2D 3\z 4\T 5\alpstt
    switch sweepsele
        case 1
            input_=[5:1:12];
        case 2
            input_=[1:1:9];
        case 3
            input_=[5:2:19];
        case 4
            input_=[110:10:150];
            %input_=[200:-10:10];
        case 5
            input_=[1:1:10];
    end
    sztot=size(input_,2);
    for ct1=1:sztot
        close all;clc
        datname=sprintf('final_%d.mat',input_(ct1));
        input_(ct1)
        figure;
        hold on
        if (0)
            if machineselc
                cd('C:\Users\Brian\Desktop\tmp\ferri')
            else
                cd('C:\Users\Brian\Desktop\tmp\ferri\m-h\neg')
            end
            load(datname);
            plot(Hextz_,mREHext_(:,3,end),'-r*')
            clear Hextz_ mREHext_
        end
        if machineselc
            cd('C:\Users\Brian\Desktop\tmp\ferri\sweepT\pos')
        else
            cd('C:\Users\a0132576\Desktop\test\ferrio\sweepT\neg')
        end
        load(datname);
        if sweepsele==5
            plot(Jc_,mREstt_(:,3,end),'-b*')
        else
            plot(Hextz_,mREHext_(:,3,end),'-b*')
        end
        %legend('init-1','init+1')
        xlabel('H_{ext}(T)');ylabel('m_{RE}')
    end
end


if(1)%for debug
    %save('tmp.mat')
    figure;
    hold on
    plot(tt_*1e12,mTM_(:,3),'linewidth',3)
    %clear mRE_
    load('final2.mat')
    plot(tt_*1e12,mTM_(:,3),'linewidth',2)
    legend('now','init')
end
if (0)%for test
    plot(tt_*1e12,reshape(mREHext_(28,3,:),[1,size(mREHext_(26,3,:),3)]),'linewidth',3)
    plot(Hextz_,mREHext_(:,3,end),'-r*')
    
    figure;
    hold on
    plot(Hextz_,mREHext_(:,3,end),'-r*')
    plot(Hextz_,mREHext_(:,3,end),'-b*')
    legend('mRE-init(-1)','mRE-init(+1)')
    title('T=150')
    xlabel('H_{ext}(K)');ylabel('mRE(z)')
end
if (0) %calc FM, e.g.TM
    if (1)%compare with experimental data[1]
        hold on;
        plot(T_,mTMFM_,'-b*');
        load('table3.mat')
        plot(t3fekx,t3feky/(t3feky(1)),'-o');
        legend('simulation','experiment')
    elseif (0)% compare with different condition.
        %save('diff-x0.mat','T_','mTM_')
        hold on
        plot(T_,mTMFM_,'b*','linewidth',3);
        clear T_ mTM_
        load('diff-x0.mat')
        plot(T_,mTM_,'r*','linewidth',2);
        legend('now','ori')
        title('FM-TM')
    else
        plot(T_,mTMFM_,'b*','linewidth',3);
    end
    xlabel('T(K)');ylabel('m')
else%RE
    if (1)%compare with experimental data[3]
        hold on;
        plot(T_,mREFM_,'b*');
        load('fig13.mat')
        plot(fig13x,fig13y/fig13y(1),'-o');
        legend('simulation','experiment')
    elseif (0)% compare with different condition.
        %save('diff-x1.mat','T_','mRE_')
        hold on
        plot(T_,mREFM_,'b*','linewidth',3);
        clear T_ mTM_
        load('diff-x1.mat')
        plot(T_,mRE_,'r*','linewidth',2);
        legend('now','ori')
        title('FM-RE')
    else
        plot(T_,mREFM_,'b*','linewidth',3);
        
    end
    xlabel('T(K)');ylabel('m')
end



