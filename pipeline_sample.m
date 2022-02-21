function millimetro_sampleCodeMultiple()
c0 = 1/sqrt(8.85e-12*4.*pi.*1e-7);
%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
CurPath = pwd();

% Brd         =   TinyRad(); %Radar parameters from Analog Devices
minRD=1000;
maxRD=-1000;
minR=1000;
maxR=-1000;
minAng=1000;
maxAng=-1000;

bkndSubtract=0;

filePre='./Data/';
c=3e8;
for experiments=1:1%5
    if experiments==1 %Set of mobility tests
        arucoGT=1;
        arucoInd=1;
        maxRang=50;
        foldern='2020-8-15-mobility';
        fileName=dir([filePre foldern '/*.mat']);
        [GT,GTAoA,xGT,yGT,GT_time]=aruco_calc([filePre foldern '/']);
        
        startFrame=1;
        frameNumOrg=1;
        switchPer=[3500]*10^-6; %625, 500
        modDuty=50;
        vel=[0];
    elseif experiments==2 %Set of tests showing angle movement, no gt
        arucoGT=0;
        maxRang=8;
        arucoInd=0
        foldern='2021-angle-tests';
        fileName=dir([filePre foldern '/*.mat']);
        GT=zeros(23,1);
        GTAoA=zeros(23,1);
        [xGT,yGT]=pol2cart(deg2rad(GTAoA+90),GT');
        xGT=xGT*-1;
        
        startFrame=1;
        frameNumOrg=1;
        switchPer=[625]*10^-6;
        modDuty=50;
        vel=[0];
    elseif experiments==3 %Set of static tests 
        arucoGT=1;
        maxRang=8;
        arucoInd=1;
        foldern='2020-8-14-indoorLoc';
        fileName=dir([filePre foldern '/*.mat']);
        GT=zeros(31,1);
        GTAoA=zeros(31,1);
        [xGT,yGT]=pol2cart(deg2rad(GTAoA+90),GT');
        xGT=xGT*-1;
        
        startFrame=1;
        frameNumOrg=1;
        switchPer=[1500]*10^-6; %625, 500
        modDuty=50;
        vel=[0];
    elseif experiments==4 %Tags within the same range bin, no gt
        arucoGT=0;
        maxRang=8;
        arucoInd=0
        foldern='same-range';
        fileName=dir([filePre foldern '/*.mat']);
        GT=zeros(23,1);
        GTAoA=zeros(23,1);
        [xGT,yGT]=pol2cart(deg2rad(GTAoA+90),GT');
        xGT=xGT*-1;
        
        startFrame=1;
        frameNumOrg=1;
        switchPer=[750 500]*10^-6; %625, 500 %750 500
        modDuty=50;
        vel=[0];
    end
    
    tagRangeErr=[];tagAoAErr=[];rangeAll=[];angleAll=[];xAll=[];yAll=[];
    
    fileName=dir([filePre foldern '/*.mat']);
    warning('off', 'MATLAB:colon:nonIntegerIndex');
    
    switchPer=unique(switchPer);
    tagHeatmap=[];
    for jj=5:size(fileName,1) %Change file index 
        disp("Running new file: ");
        disp(jj)
        %% Radar parameter retrieval and creation (Keep collapsed unless modifying)
        if(1) 
            fileInd=jj;
            fileName(fileInd).name
            RadData=load([filePre foldern '/' fileName(fileInd).name]);
            DataAll=RadData.Data;
            RadData.dtime.TimeZone='America/New_York';
            numFrame=size(DataAll,2);%For jj%size(DataAll,2); %TODO: set the number of frames to be executed (all or a handful of frames, e.g. 10)
            frameTime=RadData.dtime;

            Cfg=RadData.Cfg;
            
            if(length(Cfg.Seq) == 1)
                num_tx = 1;
            else
                num_tx = 2;
            end
            
            N=RadData.N; %number of samples
            fs=RadData.fs; %sampling freq
            num_rx = RadData.NrChn;
            NrChn = num_rx * num_tx; %number of receivers
            CalDat=RadData.CalDat;%ones(4,1);%
            measRounds=RadData.measRounds;
            DataAll(1:N:(size(DataAll,1)), :, :)=[]; %remove the chirp number
            antenna_dist=0.006217;%0.00612;

            %sincFunBased
            modF=1./(2*[switchPer]);

            % Processing of range profile
            Win2D           =   repmat(hanning(N-1),1,Cfg.FrmMeasSiz,NrChn);
            ScaWin          =   sum(Win2D(:,1,1));
            NFFT            =   2^10;%RadData.N-1;% 10

            NFFTVel         =   2^8;%size(RadData.Data,1)/RadData.N;%
            

            kf              =   (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp;%Cfg.Perd;%
            fc              =   (Cfg.fStop + Cfg.fStrt)/2;
            lambda          =   c/fc;                
            vRange          =   [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf);
            
            RMin            =   0.1;
            RMax            =   min(maxRang,((Cfg.N/Cfg.TRampUp)*c0)/(4*(250e6/Cfg.TRampUp)));

            [Val RMinIdx]   =   min(abs(vRange - RMin));
            [Val RMaxIdx]   =   min(abs(vRange - RMax));
            vRangeExt       =   vRange(RMinIdx:RMaxIdx);
            rOffset         =   length(vRangeExt(vRangeExt<RMin+0.2)); %Range offset

            doppSize        =   Cfg.FrmMeasSiz;%min(128,Cfg.FrmMeasSiz);
            WinVel          =   hanning(doppSize);
            ScaWinVel       =   sum(WinVel);
            WinVel2D        =   repmat(WinVel.',numel(vRangeExt),1);

            vFreqVel        =   [-NFFTVel./2:NFFTVel./2-1].'./NFFTVel.*(1/Cfg.Perd);
            %vVel            =   vFreqVel*c0/(2.*fc);

            % Window function for receive channelsf
            NFFTAnt         =   255;
            WinAnt          =   hanning(NrChn);
            ScaWinAnt       =   sum(WinAnt);
            WinAnt2D        =   permute(repmat(WinAnt,1,Cfg.FrmMeasSiz,numel(vRangeExt)),[3,2,1]);
            %WinAnt2D        =   repmat(WinAnt.',numel(vRangeExt)*Cfg.FrmMeasSiz,1);
            vAngDeg         =   asin(2*[-NFFTAnt./2:NFFTAnt./2-1].'./NFFTAnt)./pi*180;


            % Calibration data
            mCalData        =   permute(repmat(CalDat(1:NrChn),1,Cfg.FrmMeasSiz,N-1),[3 2 1]);

            % Positions for MeasIdxpolar plot of cfost function
            vU              =   linspace(-1,1,NFFTAnt);
            [mRange , mU]   =   ndgrid(vRangeExt,vU);
            mX              =   mRange.*mU;
            mY              =   mRange.*cos(asin(mU));
        end
        
        
       %% Initial Processing and FFTs
        disp("Beginning initial FFTs...");
        RPExtAll=zeros(numel(vRangeExt),Cfg.FrmMeasSiz,numFrame,NrChn);
        RDAll=zeros(numel(vRangeExt),NFFTVel,numFrame,NrChn);
        
        for MeasIdx = 1:1:numFrame-startFrame+1%startFrame
            Data        =  reshape(squeeze(DataAll(:,MeasIdx,:)),N-1,[], NrChn);
                        
            if(bkndSubtract) %Keep off
                Data=Data-Data(:,1,:);
            end
            
            % Calculate range profile including calibration 
            RP          =   2*fft(Data.*Win2D.*mCalData,NFFT,1).*(7.5989e-6)/ScaWin; %Replaced Brd.FuSca with constant
            RPExt       =   RP(RMinIdx:RMaxIdx,:,:);

            
            RD          =   fft(RPExt.*WinVel2D, NFFTVel, 2)./ScaWinVel;
            RDA=RD;%cat(2,RDA,RD);
            
            RPExtAll(:,:,MeasIdx-startFrame+1,:)=RPExt;
            RDAll(:,:,MeasIdx-startFrame+1,:)=RDA;
        end
        disp("Done. ");
       %% Find Range Bins
        disp("Beginning range bin search...");
        for MeasIdx = 1:1:numFrame-startFrame+1
            
            RDA   = squeeze(RDAll(:, :, MeasIdx-startFrame+1, :));

            % define the signal for correlation matching
            sig=squeeze(abs(RDA(:,:,:)));%squeeze(max(RDA(:,:,:),[],3)); %get the max value across 4 channels
            
            % normalize the signal
            normA = sig - min(sig(:));
            normSig = normA ./ max(normA(:));
            normSig = sum(normSig, 3);
            
            %Calculate background noise across frequency bins
            doppNoise=sum(sum(normSig(:,:,:),3),1);
            convW=17;
            doppNoise=conv(doppNoise, triang(convW)); %smoothing
            doppNoise=doppNoise(ceil(convW/2):end-floor(convW/2));
            doppNoise=normalize(doppNoise,2,'range');
            
            %static uses values 10 and -0.2
            %mobi using values 400*(-0.03)
            sigmoid(MeasIdx-startFrame+1, :) = 1./(1+exp(400*(-0.03+doppNoise)));
            
            %vel=[-1:0.05:1]; %for static cases
            vel=[-7.2:0.05:7.2]; %for mobility (more processing time)
            normsumR2d=zeros(length(modF), length(vel), length(vRangeExt));
            
            %Sweep each filter across results
            for vind=1:length(vel)
                for k=1:length(modF)
                    expDopFinal=calculateFilter_wMobility(Cfg,modF(k),modDuty,NFFTVel,vRangeExt, fc, c,fs,vel(vind));
                    expDopFinal(isnan(expDopFinal)) = 0;
                    expDopFinal = expDopFinal ./ (sum(expDopFinal(1, :)) + 0.0001) *3; %Want to bias away from a matched filter with more peaks/just larger
                    expDopSigmoid=sigmoid(MeasIdx-startFrame+1, :).*expDopFinal;
                    normexpDopFinal=expDopSigmoid;
                    
                    normR2=squeeze(sum(abs(normSig).*repmat(abs(normexpDopFinal),1,1,size(normSig,3)),2));
                    normsumR2d(k, vind, :)=normR2;
                end
            end
            
            newFinal=zeros(length(modF), 1);
            velocities=zeros(length(modF), 1);
            newSelR=zeros(length(modF), 1);
            
            %Find maximum value and record
            for tagNum=1:length(modF)
                [peakInds, newFinal(tagNum) ]=findPeaks2D(squeeze(normsumR2d(tagNum,:,rOffset:end)), 1);
                velocities(tagNum)=vel(peakInds(1));
                newSelR(tagNum)=peakInds(2) + rOffset;
            end
            
            finalVal=newFinal;
            selR=newSelR;
            tagLocationCorr(MeasIdx-startFrame+1, :)= squeeze(normsumR2d(1, find(vel==0), selR-2:selR+2));
            
            selRAll(jj, MeasIdx-startFrame+1, :) = newSelR;
            rangeAll(jj,MeasIdx-startFrame+1, :)=vRangeExt(newSelR-rOffset);
            velocityAll(jj,MeasIdx-startFrame+1, :)=velocities;
            corrAll(jj,MeasIdx-startFrame+1, :)=newFinal;
        end
        
        selROverall(jj, :) = round(mean(selRAll(jj, :, :), 2));
        rangeOverall(jj, :) = mean(rangeAll(jj, :, :), 2);
        velocityOverall(jj, :) = mean(velocityAll(jj, :, :), 2);
        disp("Done. ");
       %% Find Angle from Radar
        disp("Beginning AoA search... ");
        for tagNum = 1:length(modF)

            for MeasIdx = 1:1:numFrame-startFrame+1
                
                %Recreate filter
                expDopFinal = calculateFilter_wMobility(Cfg, modF(tagNum), modDuty, NFFTVel, vRangeExt, fc, c, fs, velocities(tagNum));
                
                %retrieve fft across antennas
                RDA=squeeze(RDAll(selRAll(jj,MeasIdx-startFrame+1,tagNum),:, MeasIdx-startFrame+1, :));
                
                %take range bin with tag
                RDAngle = fftshift(fft(RDA.*squeeze(WinAnt2D(1,1,:))', NFFTAnt, 2)/ScaWinAnt,2);
                rangeBin(:,:)=squeeze(abs(RDAngle));
                
                %apply sigmoid
                normexpDopFinal=expDopFinal(1,:).*sigmoid(MeasIdx-startFrame+1, :);
                    
                %Find maximum correlation
                normR2angle=squeeze(sum(rangeBin.*normexpDopFinal', 1));
                [~, aoa_ind]=max(normR2angle);
                
                %find correct angle bin
                aoa_f=vAngDeg(aoa_ind);

                angleAll(jj,MeasIdx-startFrame+1,tagNum)=aoa_f;
                angleIndAll(jj,MeasIdx-startFrame+1,tagNum)=aoa_ind;
            end
        end

        angleOverall(jj, :) = mean(angleAll(jj, :, :), 2);
        selAngleOverall(jj, :) = round(mean(angleIndAll(jj, :, :), 2));
        disp("Done. ");
        
        
       %% Do GT Calculations
       disp("Beginning GT Calcualtions...");
       for tagNum = 1:length(modF)
           for MeasIdx = 1:1:numFrame-startFrame+1
               if(arucoGT)
                    idx=findGT(frameTime(MeasIdx),GT_time);
                    if idx~=-1                       
                        GT_distAruco(arucoInd)=GT(idx);
                        GT_AoAAruco(arucoInd)=GTAoA(idx);
                        GT_xAruco(arucoInd)=xGT(idx);
                        GT_yAruco(arucoInd)=yGT(idx);
                        GT_timeAruco(arucoInd)=GT_time(idx);
                        tagRangeErr(jj,MeasIdx-startFrame+1, tagNum)=abs(rangeAll(jj, MeasIdx, tagNum) - GT(idx));
                        tagAoAErr(jj,MeasIdx-startFrame+1, tagNum)=abs(angleAll(jj, MeasIdx, tagNum) - GTAoA(idx));
                        arucoInd=arucoInd+1;
                    end
                end
           end
       end
       disp("Done. ");
    end
    end
end

%% Filter Function Generator
function expDopFinal=calculateFilter_wMobility(Cfg,modF,modDuty,NFFTVel,vRangeExt,fc,c,fs,vel)
t=[0:Cfg.Perd:Cfg.FrmMeasSiz*Cfg.Perd];%oversampling to get rid of the noise
expDopFinal=zeros(length(vRangeExt),NFFTVel);

doppShift=fft(cos(2*pi*(2*vel*fc/c).*t), NFFTVel,2);
doppShift=doppShift(1:NFFTVel/2) ;
doppShift = doppShift/max(doppShift);
doppShift(doppShift < 0.6) = 0;

undopp_square=fft(square(2*pi*(modF.').*t(1:end-1),modDuty),NFFTVel,2);
expdopt=abs(undopp_square);
expdopt=expdopt./max(expdopt);
expdopt(expdopt<0.1)=0;
expdopt=expdopt./sum(expdopt)*100; %Normalize such that each filter will have the same 'amount' of autocorrelation

if(vel>0)
    doppShiftPos=[doppShift zeros(1,NFFTVel/2)];
    sq_wav=conv(doppShiftPos, expdopt);
    expdopt=abs(sq_wav(1:NFFTVel));
elseif(vel<0)
    doppShiftNeg=flip([doppShift  zeros(1,NFFTVel/2)]);
    sq_wav=conv(doppShiftNeg, expdopt);
    expdopt=abs(sq_wav(NFFTVel:end));
else
    expdopt=abs(expdopt);
end

expdopt=expdopt./max(expdopt);
expdopt(expdopt<0.1)=0;

expDopFinal=repmat(expdopt, length(vRangeExt), 1);
end

function idx=findGT(frameTime,GT_time)
[val,idx]=min(abs(frameTime-GT_time));
if val>=seconds(1)
    idx=-1;
end
end