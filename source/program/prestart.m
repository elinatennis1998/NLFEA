% Tim Truster
% 07/11/2013

% Load algorithm library data
palgolib

%% Ask for restart file name and load the file
if strcmp(batchinter,'batch') % batch mode
    
    ibatchname = batchname;
    
    if exist('rname','var')
        clear load
        load(rname);
    else
        error('rname file does not exist')
    end
    batchinter = 'batch';
    NCR = 3;
    
    % Run batchname file driver again to reload any overwritten values
    run(ibatchname)
    
else
    
    % Store previous file names for restart file driver
    ipathname = pathname;
    ifilename = filename;
    
    if ~exist('rname','var')
    rname = uigetfile('*.m','Select the NLFEA restart file');
    end
    if isempty(rname)
        error('Must select a restart file')
    end
    clear load
    load(rname)
    batchinter = 'inter';
    NCR = 3;
    
    % Run restart file driver again to reload any overwritten values
    run([ipathname ifilename(1:end-2)])

end

if currstep == 'y'
            
elseif currstep == 'n'
            
    % Load data from desired step; note that history variable
    % storage is required for this functionality.
    
    if stepin ~= step

    step = stepin;
    step = step-1; %correction on 8/21/13 to data storage

    % Load data
    lamda = TimeList(step+1);
    if DHist == 1
       ModelDx(NDOFT2(DOFa))  = DispList(DOFa+step*ndf*numnp);
       gBC(NDOFT2(DOFi)-neq)  = DispList(DOFi+step*ndf*numnp);
       ModelDxn_1(NDOFT2(DOFa))  = ModelDx;
       gBC_n(NDOFT2(DOFi)-neq)  = gBC;
    else
        error('Displacement history is unavailable')
    end
    if ismember(transient,transalgos)
    if VHist == 1
       ModelVx(NDOFT2(DOFa))  = VeloList(DOFa+step*ndf*numnp);
       ModelVxn_1(NDOFT2(DOFa))  = ModelVx;
    else
        error('Velocty history is unavailable')
    end
    if AHist == 1
       ModelAx(NDOFT2(DOFa))  = AcceList(DOFa+step*ndf*numnp);
       ModelAxn_1(NDOFT2(DOFa))  = ModelAx;
    else
        error('Acceleration history is unavailable')
    end
    end

    if ~exist('hisover','var')
        hisover = 0;
    end
    if hrstore
        hrvec = hrmat(:,step+1);
    else
        if hisover
            disp('Warning: history data cannot be loaded') %#ok<UNRCH>
        else
            error('History data DNE; if history is not needed, set hisover=1')
        end
    end

    step = step+1; %correction on 8/21/13 to data storage
    
    end
        
end