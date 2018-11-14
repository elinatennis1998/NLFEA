

%Tim Truster
%Non-Linear Finite Element Program, ver2
%Created: 04/2013

format compact

if (~exist('batchinter','var'))
    batchinter = 'inter';
end

if (~exist('NCR','var'))
    NCR = 1; % 1=New, 2=Continuation, 3=Restart
end

%% Launch Program: Select Either new analysis, Continuation, or Restart

if strcmp(batchinter,'batch') % batch mode
    
    
    if NCR == 1 % New simulation

        if exist('batchname','var')
            if strcmp(batchname(end-1:end),'.m')
                batchname = batchname(1:end-2);
            end
            run(batchname);
        else
            error('batchname file does not exist')
        end

        if NCR ~= 1
            error('NCR~=1 after reading file; this file may not be an input file')
        end
        batchinter = 'batch';
        
        % Read major problem size parameters and solution options
        pstart
        
        % Initialize default parameters
        pdefault
        
    elseif NCR == 2 % Continuation

        if exist('batchname','var')
            if strcmp(batchname(end-1:end),'.m')
                batchname = batchname(1:end-2);
            end
            run(batchname);
        else
            error('batchname file does not exist')
        end

        if NCR ~= 2
            error('NCR~=2 after reading file; this file may not be a continuation file')
        end
        batchinter = 'batch';
        
        % Initialize default parameters
        pdefault
        
    elseif NCR == 3 % Restart

        if exist('batchname','var')
            if strcmp(batchname(end-1:end),'.m')
                batchname = batchname(1:end-2);
            end
            run(batchname);
        else
            error('batchname file does not exist')
        end

        if NCR ~= 3
            error('NCR~=3 after reading file; this file may not be a restart file')
        end
        batchinter = 'batch';

        % Request input for restart file and step number
        prestart
        
        % Initialize default parameters
        pdefault

    else

        error('invalid choice for NCR')
        
    end    
    
else % interactive - prompt for filename
    
    if NCR == 1 % New simulation

    %% Ask for input file name and load the file
        if exist('filename.mat','file')
        clear load
        load('filename.mat', 'defaultname','defaultpath')
        else
        defaultname = 'DNE';
        end
        fprintf('Provide name of input file (default is %s)\n',defaultname);
        [filename,pathname] = uigetfile('*.m','Select the NLFEA input file');
        if isnumeric(filename)
             pathname = defaultpath;
             filename = defaultname;
        end
        defaultname = filename;
        defaultpath = pathname;
        save('filename.mat','defaultname','defaultpath');
        run([pathname filename(1:end-2)])
        load('filename.mat') % reload filename so that the workspace knows which file it was

        if NCR ~= 1
            error('NCR~=1 after reading file; this file may not be an input file')
        end
        batchinter = 'inter';
        
        % Read major problem size parameters and solution options
        pstart
        
        % Initialize default parameters
        pdefault

    elseif NCR == 2 % Continuation

    %% Ask for continuation input file name and load the file
        [filename,pathname] = uigetfile('*.m','Select the NLFEA continuation file');
        if isnumeric(filename) % User cancelled
            error('Must select a continuation file')
        end
        run([pathname filename(1:end-2)])

        if NCR ~= 2
            error('NCR~=2 after reading file; this file may not be a continuation file')
        end
        batchinter = 'inter';
        
        % Initialize default parameters
        pdefault

    elseif NCR == 3 % Restart

    %% Ask for restart input file name and load the file
        [filename,pathname] = uigetfile('*.m','Select the NLFEA restart file');
        if isnumeric(filename) % User cancelled
            error('Must select a restart file')
        end
        run([pathname filename(1:end-2)])

        if NCR ~= 3
            error('NCR~=3 after reading file; this file may not be a restart file')
        end
        batchinter = 'inter';

        % Request input for restart file and step number
        prestart
        
        % Initialize default parameters
        pdefault

    else

        error('invalid choice for NCR')

    end
    
end

perror

tic

%% FE array initialization
if NCR == 1

    lamda = 0; %time
    step = 0;
    iter = 0;
    failFEelem = 0;
    initializeFE

    %Set up sparsity pattern, count number of lin/nonlin materials
    FormKs

    % Initialize solution vectors; compute any required initial arrays

    ModelFx = zeros(nieq,1);
    s_del_ModelDx = zeros(neq,1);
    gBC = zeros(nieq,1);
    gBC_n = gBC;
%     ModelVx = zeros(neq,1);
%     ModelAx = zeros(neq,1);
%     ModelDxn_1 = ModelDx;
%     ModelVxn_1 = ModelVx;
%     ModelAxn_1 = ModelAx;
    % dt = s_del_a;
    initia = 1; %flag to prescribe an elastic predictor step
    NDOFT2 = NDOFT';
    DOFa = find((NDOFT2<=neq)&(NDOFT2>0));
    DOFi = find(NDOFT2>neq);

    isw = 40;
    FormFE

    if Compt == 1
    %     PatchTimes
    %     IterTimes
    %     StepTimes
    %     SolveTimes
    end
    %History flags for various output arrays
    TimeList = zeros(1,datastep);
    if IHist == 1
        IterList = zeros(3+itermax,datastep);
    end
    if DHist == 1
        DispList = zeros(ndf,numnp,datastep);
        DispList0 = zeros(ndf,numnp);
        DispList0(DOFa) = ModelDx(NDOFT2(DOFa));
    end
    if VHist == 1
        VeloList = zeros(ndf,numnp,datastep);
        VeloList0 = zeros(ndf,numnp);
        VeloList0(DOFa) = ModelVx(NDOFT2(DOFa));
    end
    if AHist == 1
        AcceList = zeros(ndf,numnp,datastep);
        AcceList0 = zeros(ndf,numnp);
        AcceList0(DOFa) = ModelAxn_1(NDOFT2(DOFa));
    end


    if FHist == 1
        ForcList = zeros(ndf,numnp,datastep);
    end


    % Initialize lists for output data
    poutinit
    if LieHist == 1
        LieHList = zeros(numLie,numnp,nummat,datastep);
        numSubnormvals = 6;
        SubnormsT = zeros(numSubnormvals,itermax,datastep); % Initial value for norms
        SubnormsL = zeros(numSubnormvals,itermax*(numsubcycles+1),datastep); % norm of changes to nodal values of Re
%         SubStreL = zeros(100,stepmax); % norms of changes to Cauchy stress
%         SubTauL = zeros(100,stepmax); % norms of changes to tau_tilde, the hardening variable
    end

elseif NCR == 2 || NCR == 3
    
    % Add more external loads
    if (exist('newloads','var') && newloads)
        ploadi_new
    end
    
    if transient == 0 && Fcount == 0 && norm(Fc1+Fc1np) == 0
        b = -1;
    end

end


%% Problem Initialization
if NCR == 1 || NCR == 2
switch transient
    
    case -1 %linear static analysis
        
        step = 1; lamda = 1;
        GetAllLoads
        GetAllBCs
        gBC_n = zeros(nieq,1);
        
        %Assemble Stiffness Routine
        initializeLinKF = 1;
        isw = 3;
        FormFE
        initializeLinKF = 0;
    
    case 0 %initial initialization quasi-static 
        
        lamdamax = max(mults)+1;
        
        if NCR == 1
        
        lamda = 0;
        GetAllLoads
        gBC = zeros(nieq,1);
        gBC_n = gBC;
        
        %Assemble Stiffness Routine
        initializeLinKF = 1;
        isw = 3;
        FormFE
        initializeLinKF = 0;
        
        Fdisp = -Kdf1*ModelDc;
        
        end
        
        if printRnorm
            fprintf('Rnorm          | log_10(Rnorm)   | log_10(Rnorm) - 2*log_10(Rnorm1)\n')
        end
        
    case 1 %Alpha-Newmark Method d-form
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        F1n_1 = Nalpha*Fext1;
        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        coeffk = 1;
        coeffm = 1;
        
        %Compute M
        % Mass lumping should not be used for total energy conservation
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        initializeLinKF = 1;
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
%         Mdd22 = diag(sum(Mdd11))*1.5;
        ModelAxn_1 = Mdd11\Fd1;
%         ModelAxn_1 = Mdd22\Fd1;
    if AHist == 1
        AcceList0(DOFa) = ModelAxn_1(NDOFT2(DOFa));
    end
        
        %Compute M
        isw = 5;
        FormFE

        initializeLinKF = 1;
        isw = 21;
        FormFE
        initializeLinKF = 0;
        
        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = 1;
        coeffv = Ngamma/(Nbeta*tstep);
        coeffa = 1/(Nbeta*tstep^2);
        % Linear Matrix coefficients
        coeffkl1 = -Nalpha;
        coeffml1 = 0;
        coeffkl2 = (1 + Nalpha);
        coeffml2 = 1;
        % Matrix coefficients
        coeffk = (1 + Nalpha);
        coeffm = 1/(Nbeta*tstep^2);
        Kprime = coeffm*Mdd11;
        if numbernonlinear == 0 && exist('factorize','file')
        Kn = Kprime + coeffk*Kdd11;
%         [KnR,perr,Schol] = chol(Kn,'matrix');
        KnR = factorize(Kn);
        end
        
        if enercons == 1
            isw = 12;
            FormFE
            if numberlinear > 0
            SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
            end
            En = SysEner
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end

    case 2 %Alpha-Newmark Method a-form
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        F1n_1 = Nalpha*Fext1;
        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        coeffk = 1;
        coeffm = 1;
        
        %Compute M
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        initializeLinKF = 1;
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
        ModelAxn_1 = Mdd11\Fd1;
    if AHist == 1
        AcceList0(DOFa) = ModelAxn_1(NDOFT2(DOFa));
    end
        
        initializeLinKF = 1;
        isw = 21;
        FormFE
        initializeLinKF = 0;

        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = Nbeta*tstep^2;
        coeffv = Ngamma*tstep;
        coeffa = 1;
        % Linear Matrix coefficients
        coeffkl1 = -Nalpha;
        coeffml1 = 0;
        coeffkl2 = (1 + Nalpha);
        coeffml2 = 1;
        % Matrix coefficients
        coeffk = (1 + Nalpha)*Nbeta*tstep^2;
        coeffm = 1;
        Kprime = coeffm*Mdd11;
        if numbernonlinear == 0 && exist('factorize','file')
        Kn = Kprime + coeffk*Kdd11;
%         [KnR,perr,Schol] = chol(Kn,'matrix');
        KnR = factorize(Kn);
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end
        
    case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        F1n_1 = Nalpha*Fext1;
        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        
        %Compute M
        % Mass lumping should not be used for total energy conservation
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        initializeLinKF = 1;
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
%         Mdd22 = diag(sum(Mdd11))*1.5;
        ModelAxn_1 = Mdd11\Fd1;
%         ModelAxn_1 = Mdd22\Fd1;
        
        %Compute M
        isw = 5;
        FormFE

        %Compute Kprime, Set coefficients
        isw = 21;
        FormFE
        initializeLinKF = 0;
        
        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = 1;
        coeffv = Ngamma/(Nbeta*tstep);
        coeffa = 1/(Nbeta*tstep^2);
        % Linear Matrix coefficients
        coeffkl1 = -Nalpha;
        coeffml1 = 0;
        coeffkl2 = (1 + Nalpha);
        coeffml2 = 1;
        % Matrix coefficients
        coeffk = (1 + Nalpha);
        coeffm = 1/(Nbeta*tstep^2);
        Kprime = coeffm*Mdd11;
        
        if enercons == 1
            isw = 12;
            FormFE
            if numberlinear > 0
            SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
            end
            En = SysEner
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end
        
        if exist('lenerHist','var') && lenerHist == 1
            lenerList = zeros(datastep,1);
        end
        
    case 4 %Energy Conservation by LM according to Hughes, Trans of ASME
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        F1n_1 = Nalpha*Fext1;
        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        
        %Compute M
        % Mass lumping should not be used for total energy conservation
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
%         Mdd22 = diag(sum(Mdd11))*1.5;
        ModelAxn_1 = Mdd11\Fd1;
%         ModelAxn_1 = Mdd22\Fd1;
        
        %Compute M
        isw = 5;
        FormFE

        %Compute Kprime, Set coefficients
        initializeLinKF = 1;
        isw = 21;
        FormFE
        initializeLinKF = 0;
        coeff0 = 2/(tstep^2);
        coeff6 = tstep;
        coeff7 = 1/2*(1 - 2*beta)*tstep^2;
        coeff8 = (1 - gamma)*tstep;
        Kprime = coeff0*Mdd11;
        
        if enercons == 1
            isw = 12;
            FormFE
            SysEner = SysEner + 1/2*ModelVx'*Mdd11*ModelVx;
            En = SysEner
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end
        
    case 5 % Generalized Alpha-Newmark Method d-form
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        Fd1 = Fext1;
        F1n_1 = 0*Fext1;
        F1n = F1n_1;
%         % FEAP does not initialize the accel for HHT method
%         coeff0 = 1;
%         coeff1 = 1;
%         coeff2 = 1;
%         coeff3 = 1;
%         coeff4 = 1;
%         coeff5 = 1;
%         coeff6 = 0;
%         coeff7 = 0;
%         coeff8 = 0;
%         
%         %Compute M
% %         lumping = 1; % Add these two lines to do an initial lump like FEAP
%         isw = 5;
%         FormFE
% %         lumping = 0;
% 
%         %Compute Residual
%         isw = 6;
%         FormFE
%         
%         %Compute initial accelerations
%         ModelAxn_1 = Mdd11\Fd1;
%         % FEAP does not initialize the accel for HHT method
        
        %Compute M
        initializeLinKF = 1;
        isw = 5;
        FormFE

        %Compute Kprime, Set coefficients
        isw = 21;
        FormFE
        initializeLinKF = 0;
        
        % Predictors
        coeffdd = 1;
        coeffdv = 0;
        coeffda = 0;
        coeffvd = 0;
        coeffvv = 1 - Ngamma/Nbeta;
        coeffva = tstep*(1 - Ngamma/(2*Nbeta));
        coeffad = 0;
        coeffav = -1/(Nbeta*tstep);
        coeffaa = (1 - 1/(2*Nbeta));
%         coeffdd = 1;
%         coeffdv = tstep;
%         coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
%         coeffvd = 0;
%         coeffvv = 1;
%         coeffva = (1 - Ngamma)*tstep;
%         coeffad = 0;
%         coeffav = 0;
%         coeffaa = 0;
        % Correctors
        coeffd = 1;
        coeffv = Ngamma/(Nbeta*tstep);
        coeffa = 1/(Nbeta*tstep^2);
        % Linear Matrix coefficients
        coeffkl1 = Nalphaf;
        coeffml1 = Nalpham;
        coeffkl2 = (1 - Nalphaf);
        coeffml2 = (1 - Nalpham);
        % Matrix coefficients
        coeffk = (1 - Nalphaf);
        coeffm = (1 - Nalpham)/(Nbeta*tstep^2);
        Kprime = coeffm*Mdd11;
        if numbernonlinear == 0 && exist('factorize','file')
        Kn = Kprime + coeffk*Kdd11;
%         [KnR,perr,Schol] = chol(Kn,'matrix');
        KnR = factorize(Kn);
        end
        
        if enercons == 1
            isw = 12;
            FormFE
            if numberlinear > 0
            SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
            end
            En = SysEner
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end
        
    case 6 % Multiscale Dynamics
        
        % Set flag for using acceleration form (2) or displacement form (1)
        adform = 1;2;
        
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        Fd1 = Fext1;
        F1n_1 = 0*Fext1;
        F1n = F1n_1;

        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        coeffk = 1;
        coeffm = 1;
        
        %Compute M
        initializeLinKF = 1;
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        fsinit = 1;
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
        ModelAxn_1 = Mdd11\Fd1;
        
        % Initialize FS accel
        isw = 40;
        FormFE
        
        %Compute M
        fsinit = 0;
        isw = 5;
        FormFE

%         %Compute Kprime, Set coefficients % removed soas to not evolve the fine scales before energy calculation
%         initializeLinKF = 1;
%         isw = 21;
%         FormFE
        initializeLinKF = 0;
        
        if adform == 1 %d-form
        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = 1;
        coeffv = Ngamma/(Nbeta*tstep);
        coeffa = 1/(Nbeta*tstep^2);
        % Linear Matrix coefficients
        coeffkl1 = 0;
        coeffml1 = 0;
        coeffkl2 = 1;
        coeffml2 = 1;
        % Matrix coefficients
        coeffk = 1;
        coeffm = 1/(Nbeta*tstep^2);
        else %a-form
        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = Nbeta*tstep^2;
        coeffv = Ngamma*tstep;
        coeffa = 1;
        % Linear Matrix coefficients
        coeffkl1 = 0;
        coeffml1 = 0;
        coeffkl2 = 1;
        coeffml2 = 1;
        % Matrix coefficients
        coeffk = Nbeta*tstep^2;
        coeffm = 1;
        end
        
        if enercons == 1
            isw = 12;
            FormFE
            En = SysEner
            maxen = SysEner;
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end
    
    case 7 %initial initialization for arc-length
        
        if NCR == 1
        
        step = 0; lamda = 1;
        GetAllLoads
        gBC = zeros(nieq,1);
        gBC_n = gBC;
        if Fcount == 0 && norm(Fext1) == 0
            b = -1;
        end
        
        %Assemble Stiffness Routine
        initializeLinKF = 1;
        isw = 3;
        FormFE
        initializeLinKF = 0;
        
        Fdisp = -Kdf1*ModelDc;
        
        K0 = diag(Kdd11);
        if b > 0 %force control or arc-length
            clear('det')
            detKn_1 = det(Kdd11);
%             s_del_lamda_n0 = sign(detKn_1)*s_del_a; %goes negative for mixed element
            s_del_lamda_n0 = s_del_a;
            s_del_lamda_n1 = s_del_lamda_n0;
            detKn_2 = detKn_1;
            q0 = Kdd11\(Fc1+Fc1np+FcU);
            c = 0;
%             for i = 1:neq
%                 c = c + q0(i)*K0(i)*q0(i);
%             end
            c = c + sum(q0.*K0.*q0);
            c = (1-b)/c;
        else %disp control
            s_del_lamda_n1 = s_del_a;
        end
        
        end
        
        if printRnorm
            fprintf('Rnorm          | log_10(Rnorm)   | log_10(Rnorm) - 2*log_10(Rnorm1)\n')
        end
        
    case 8 % Thermodynamics
        
        % Get temp and disp dofs
        NDOFT3 = diag([0 0 1])*NDOFT2;
%         tdofs = NDOFT2(find((NDOFT3<=neq)&(NDOFT3>0)));
        tdofs = NDOFT2((NDOFT3<=neq)&(NDOFT3>0));
        ddofs = setdiff((1:neq)',tdofs);
        
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        Fd1 = Fext1;
        F1n_1 = 0*Fext1;
        F1n = F1n_1;

        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        
        %Compute M
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        initializeLinKF = 1;
        fsinit = 1;
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
        ModelAxn_1 = Mdd11\Fd1;
        ModelVxn_1(tdofs) = ModelAxn_1(tdofs);
        ModelAxn_1(tdofs) = 0;
        
        % Initialize FS accel
        isw = 40;
        FormFE
        
        %Compute M
        fsinit = 0;
        isw = 5;
        FormFE

%         %Compute Kprime, Set coefficients % removed soas to not evolve the fine scales before energy calculation
%         initializeLinKF = 1;
%         isw = 21;
%         FormFE
        initializeLinKF = 0;
        
        coeff0 = 1/(Nbeta*tstep^2);
        coeff1 = 1;
        coeff2 = Ngamma/(Nbeta*tstep);
        coeff3 = 1;
        coeff4 = coeff1*coeff2;
        coeff5 = coeff1;
        coeff6 = tstep;
        coeff7 = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeff8 = (1 - Ngamma)*tstep;
%         Kprime = (1 - Nalpham)/(Nbeta*tstep^2)*Mdd11;
%         if numbernonlinear == 0
%         Kn = Kprime + coeff5*Kdd11;
% %         [KnR,perr,Schol] = chol(Kn,'matrix');
%         KnR = factorize(Kn);
%         end
        
        if enercons == 1
            isw = 12;
            FormFE
%             SysEner = SysEner + 1/2*ModelVx'*Mdd11*ModelVx;
            En = SysEner
            maxen = SysEner;
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end

    case 10 % Explicit dynamics a-form
        
        Nbeta = 0;
        Ngamma = 1/2;
        Nalpha = 0;
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        
        %Compute M
        lumping = 1; % Explicit method uses lumped mass
        initializeLinKF = 1;
        isw = 5;
        FormFE

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
        ModelAxn_1 = Mdd11\Fd1;
        
        % Compute stiffness matrix for linear elements
        initializeLinKF = 1;
        isw = 21;
        FormFE
        initializeLinKF = 0;

        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = Nbeta*tstep^2;
        coeffv = Ngamma*tstep;
        coeffa = 1;
        
        if printRnorm
            fprintf('Rnorm\n')
        end
    
    case 11 % Crystal plasticity with averaging of Fe_inv
        
        lamdamax = max(mults)+1;
        
        if NCR == 1
        
        step = 0; lamda = 0;
        GetAllLoads
        gBC = zeros(nieq,1);
        gBC_n = gBC;
        
%         %Assemble Stiffness Routine
%         initializeLinKF = 1;
%         isw = 3;
%         FormFE
        initializeLinKF = 0;
%         
%         Fdisp = -Kdf1*ModelDc;
        
        end
        
        if printRnorm
            fprintf('Rnorm          | R_abs_max      | log_10(Rnorm)   | log_10(Rnorm) - 2*log_10(Rnorm1)\n')
        end

    case 12 % Forward transient Euler, based off of case 10
        
        Nbeta = 0;
        Ngamma = 1/2;
        Nalpha = 0;
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        
        %Compute M
        lumping = 1; % Explicit method uses lumped mass
        initializeLinKF = 1;
        isw = 5;
        FormFE

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
        ModelAxn_1 = Mdd11\Fd1;
        
        % Compute stiffness matrix for linear elements
        initializeLinKF = 1;
        isw = 21;
        FormFE
        initializeLinKF = 0;

        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = Nbeta*tstep^2;
        coeffv = Ngamma*tstep;
        coeffa = 1;
        
        if printRnorm
            fprintf('Rnorm\n')
        end
        
    case 13 % Attempt at GLS enforcing energy constraint
        s_del_lamda_n1 = tstep;
        step = 0;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        F1n_1 = Nalpha*Fext1;
        coeff0 = 1;
        coeff1 = 1;
        coeff2 = 1;
        coeff3 = 1;
        coeff4 = 1;
        coeff5 = 1;
        coeff6 = 0;
        coeff7 = 0;
        coeff8 = 0;
        
        %Compute M
        % Mass lumping should not be used for total energy conservation
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        initializeLinKF = 1;
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
%         Mdd22 = diag(sum(Mdd11))*1.5;
        ModelAxn_1 = Mdd11\Fd1;
%         ModelAxn_1 = Mdd22\Fd1;
        
        %Compute M
        isw = 5;
        FormFE

        %Compute Kprime, Set coefficients
        isw = 21;
        FormFE
        initializeLinKF = 0;
        
        % Predictors
        coeffdd = 1;
        coeffdv = tstep;
        coeffda = 1/2*(1 - 2*Nbeta)*tstep^2;
        coeffvd = 0;
        coeffvv = 1;
        coeffva = (1 - Ngamma)*tstep;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = 1;
        coeffv = Ngamma/(Nbeta*tstep);
        coeffa = 1/(Nbeta*tstep^2);
        % Linear Matrix coefficients
        coeffkl1 = -Nalpha;
        coeffml1 = 0;
        coeffkl2 = (1 + Nalpha);
        coeffml2 = 1;
        % Matrix coefficients
        coeffk = (1 + Nalpha);
        coeffm = 1/(Nbeta*tstep^2);
        Kprime = coeffm*Mdd11;
        
        if enercons == 1
            isw = 12;
            FormFE
            if numberlinear > 0
            SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
            end
            En = SysEner
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end
        
        if exist('lenerHist','var') && lenerHist == 1
            lenerList = zeros(datastep,1);
        end
        
    case 14 % Backward transient Euler, based off of case 5
        s_del_lamda_n1 = tstep;
        step = 0; lamda = 1;
        GetAllBCs
        gBC_n = gBC;
        lamda = 0;
        GetAllLoads
        Fd1 = Fext1;
        F1n_1 = 0*Fext1;
        F1n = F1n_1;
%         % FEAP does not initialize the accel for HHT method
%         coeff0 = 1;
%         coeff1 = 1;
%         coeff2 = 1;
%         coeff3 = 1;
%         coeff4 = 1;
%         coeff5 = 1;
%         coeff6 = 0;
%         coeff7 = 0;
%         coeff8 = 0;
%         
%         %Compute M
% %         lumping = 1; % Add these two lines to do an initial lump like FEAP
%         isw = 5;
%         FormFE
% %         lumping = 0;
% 
%         %Compute Residual
%         isw = 6;
%         FormFE
%         
%         %Compute initial accelerations
%         ModelAxn_1 = Mdd11\Fd1;
%         % FEAP does not initialize the accel for HHT method

        %Compute M
%         lumping = 1; % Add these two lines to do an initial lump like FEAP
        initializeLinKF = 1;
        isw = 5;
        FormFE
%         lumping = 0;

        %Compute Residual
        isw = 6;
        FormFE
        
        %Compute initial accelerations
        ModelVxn_1 = Mdd11\Fd1;
        
        
        %Compute M
        initializeLinKF = 1;
        isw = 5;
        FormFE

        %Compute Kprime, Set coefficients
        isw = 21;
        FormFE
        initializeLinKF = 0;
        
        % Predict v as 0
        % Predictors
        coeffdd = 1;
        coeffdv = (1-Nalpha)*tstep;
        coeffda = 0;
        coeffvd = 0;
        coeffvv = 0;
        coeffva = 0;
        coeffad = 0;
        coeffav = 0;
        coeffaa = 0;
        % Correctors
        coeffd = 1;
        coeffv = 1/(Nalpha*tstep);
        coeffa = 0;
%         % Predict v as v_n-1
%         % Predictors
%         coeffdd = 1;
%         coeffdv = tstep;
%         coeffda = 0;
%         coeffvd = 0;
%         coeffvv = 1;
%         coeffva = 0;
%         coeffad = 0;
%         coeffav = 0;
%         coeffaa = 0;
%         % Correctors
%         coeffd = 1;
%         coeffv = 1/(Nalpha*tstep);
%         coeffa = 0;
        % Linear Matrix coefficients
        coeffkl1 = 0;
        coeffml1 = 0;
        coeffkl2 = 1;
        coeffml2 = 1;
        % Matrix coefficients
        coeffk = 1;
        coeffm = 1/(Nalpha*tstep);
        Kprime = coeffm*Mdd11;
        if numbernonlinear == 0
        Kn = Kprime + coeffk*Kdd11;
%         [KnR,perr,Schol] = chol(Kn,'matrix');
        KnR = factorize(Kn);
        end
        
        if enercons == 1
            isw = 12;
            FormFE
            if numberlinear > 0
            SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
            end
            En = SysEner
        end

%         if damp == 1
%             Cdd11 = ac*Mdd11 + bc*Kdd11;
%             Kprime = Kprime + coeff4*Cdd11;
%         end
        
        if printRnorm
            fprintf('Rnorm\n')
        end
        
end
end


%% Loop for time-history of problem solution

if transient >= 0

    NR_Loop3
    
else
    
    %Solve Kd = F

    del_ModelDx = Kdd11\Fd1;

    ModelDx = ModelDx + del_ModelDx;
    
        if DHist == 1
%         for node = 1:numnp
%         %     Node_U_V(node, 1) = node;
%             for dir = 1:ndf
%                 gDOF = NDOFT(node,dir);
%                 if gDOF <= neq && gDOF > 0
%                     DispList(dir,node,step+1) = ModelDx(gDOF);
%                 elseif gDOF > 0
%                     DispList(dir,node,step+1) = gBC(gDOF - neq);
%                 end
%             end
%         end
        DispList(DOFa) = ModelDx(NDOFT2(DOFa));
        DispList(DOFi) = gBC(NDOFT2(DOFi)-neq);
        end
    if FHist == 1
        Fext1 = zeros(neq,1);
        isw = 6;
        FormFE
        
%         for node = 1:numnp
%         %     Node_U_V(node, 1) = node;
%             for dir = 1:ndf
%                 gDOF = NDOFT(node,dir);
%                 if gDOF <= neq && gDOF > 0
%                     ForcList(dir,node,step+1) = Fd1(gDOF);
%                 elseif gDOF > 0
%                     ForcList(dir,node,step+1) = Fd3(gDOF - neq);
%                 end
%             end
%         end
        ForcList(DOFa) = Fd1(NDOFT2(DOFa));
        ForcList(DOFi) = Fd3(NDOFT2(DOFi)-neq);
    end
    if SHist == 1

        StreList2 = zeros(numnp,npstr);
        Eareas = zeros(numnp,1);

        isw = 25;
        FormS2

        for stres = 1:npstr
            StreList2(:,stres) = StreList2(:,stres)./Eareas;
        end
        StreList(1:npstr,:,1) = StreList2';

    end
     
    if SEHist == 1

        StreList2 = zeros(numel,nestr);

        isw = 26;
        FormFE

        StreListE(1:nestr,:,1) = StreList2';

    end

    if expliciterr == 1

        Explicit_Error_Estimation

    end

    if impliciterr == 1

        Implicit_Error_Estimation

    end
    
end


%% Perform final data processing

Node_U_V = zeros(numnp,ndf)';

Node_U_V(DOFa) = ModelDx(NDOFT2(DOFa));
Node_U_V(DOFi) = gBC(NDOFT2(DOFi)-neq);
Node_U_V = Node_U_V';
% for node = 1:numnp
%     for dir = 1:ndf
%         gDOF = NDOFT(node, dir);
%         if gDOF > 0
%         if gDOF <= neq
%             Node_U_V(node, dir) = ModelDx(gDOF,1);
%         else
%             Node_U_V(node, dir) = gBC(gDOF - neq);
%         end
%         end
%     end
% end
