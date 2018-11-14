function AllOutput = mm10_plugin_lib(task,c_props,AllInput)
% Tim Truster
% 04/28/2016
%
% Making a task brancher for CP models, to simplify organization
% Each plugin now becomes a hardening model

    plugin = c_props.plugin;
    switch plugin
        case 1 % Ma-Roters-Raabe model
            AllOutput = mm10_plugin_mrr(task,c_props,AllInput);
        case 2 % SUVIC model
            AllOutput = mm10_plugin_suvic(task,c_props,AllInput);
        otherwise
            plugin
            error('model not implemented for this plugin')
    end

end

function AllOutput = mm10_plugin_template(task,props,AllInput)
% Plugin function for Ma-Roters-Raabe model
    switch task
        case 1 % setup_mm10_rknstr
            % AllInput = 0
            % AllOutput = num_hard
        case 2 % mm10_init_cc_hist0
            % AllInput = props, history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard),...
            %    history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard)
            % AllOutput = props, history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard), ...
            %    history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard)
        case 3 % mm10_setup
            % AllInput = props, np1, n
            % AllOutput = props, np1, n
        case 4 % mm10_formR2
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = h
        case 5 % mm10_formR2i
            error('not a function in matlab')
        case 6 % mm10_form_dbarp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = dbar
        case 7 % mm10_form_wbarp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = wbar
        case 8 % mm10_form_wp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = w
        case 9 % mm10_form_dbarpi
            error('not a function in matlab')
        case 10 % mm10_form_wpi
            error('not a function in matlab')
        case 11 % mm10_formJ11
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = dgammadtau
        case 12 % mm10_formJ12
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = dgammadtt
        case 13 % mm10_formJ21
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = estr
        case 14 % mm10_formJ22
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = etau
        case 15 % mm10_formvecs
            % AllInput = props, np1, n, stress, tt, vec1, vec2
            % AllOutput = props, np1, n, stress, tt, vec1, vec2
        case 16 % mm10_formarrs
            % AllInput = props, np1, n, stress, tt, vec1, vec2, arr1, arr2,both
            % AllOutput = arr1, arr2
        case 17 % mm10_formvecsi
            error('not a function in matlab')
        case 18 % mm10_output
            % AllInput = props, np1, n, vec1, vec2, np1.stress, np1.tau_tilde
            % AllOutput = np1.slip_incs(1:props.nslip)
        case 19 % mm10_tangent->mm10_ed_mrr
            % AllInput = props, np1, n, np1.stress, np1.tau_tilde
            % AllOutput = ed
        case 20 % mm10_tangent->mm10_dgdd_mrr
            % AllInput = props, np1, n, np1.stress, np1.tau_tilde, np1.D
            % AllOutput = dgammadd
        otherwise
            task
            error('task not implemented')
    end
end

function AllOutput = mm10_plugin_mrr(task,props,AllInput)
% Plugin function for Ma-Roters-Raabe model
    switch task
        case 1 % setup_mm10_rknstr
            % AllInput = 0
            % AllOutput = num_hard
            AllOutput = props.nslip;
        case 2 % mm10_init_cc_hist0
            % AllInput = props, history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard),...
            %    history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard)
            % AllOutput = props, history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard), ...
            %    history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard)
            [Output1, Output2, ...
                Output3] =...
                 mm10_init_mrr(AllInput{1}, AllInput{2},AllInput{3});
            AllOutput = {Output1 Output2 Output3};
        case 3 % mm10_setup
            % AllInput = props, np1, n
            % AllOutput = props, np1, n
            [Output1, Output2, ...
                Output3] = mm10_setup_mrr(AllInput{1}, AllInput{2},AllInput{3});
            AllOutput = {Output1 Output2 Output3};
        case 4 % mm10_formR2
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = h
            [~, ~, ~, ~, ~, ~, ~, h] = mm10_h_mrr(AllInput{1}, AllInput{2},...
            AllInput{3},AllInput{4},AllInput{5},AllInput{6},AllInput{7});
            AllOutput = h;
        case 5 % mm10_formR2i
            error('not a function in matlab')
        case 6 % mm10_form_dbarp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = dbar
            np1 = AllInput{2};
            vec1 = AllInput{4};
            dbar = np1.ms*vec1(1:props.nslip);
            AllOutput = dbar;
        case 7 % mm10_form_wbarp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = wbar
            np1 = AllInput{2};
            vec1 = AllInput{4};
            wbar = np1.qs*vec1(1:props.nslip);
            AllOutput = wbar;
        case 8 % mm10_form_wp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = w
%         w2 = np1.qc*vec1(1:props.nslip);
            w = zeros(3,1);
            np1 = AllInput{2};
            gam = AllInput{4}*0;
            for i = 1:props.nslip
                gam(i) = mm10_slipinc_mrr(AllInput{1}, np1, AllInput{3}, AllInput{6}, AllInput{7}, i);
                w = w + gam(i)*...
                      np1.qc(1:3,i);
            end
            AllOutput = w;
        case 9 % mm10_form_dbarpi
            error('not a function in matlab')
        case 10 % mm10_form_wpi
            error('not a function in matlab')
        case 11 % mm10_formJ11
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = dgammadtau
            arr1 = AllInput{6};
            dgammadtau = arr1(1:props.num_hard,1)';
            AllOutput = dgammadtau;
        case 12 % mm10_formJ12
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = dgammadtt
            arr2 = AllInput{7};
            dgammadtt = arr2(1:props.nslip,1:props.num_hard);
            AllOutput = dgammadtt;
        case 13 % mm10_formJ21
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = estr
            [~, ~, ~, ~, ~, ~, ~, ~, ~, estr] = mm10_estress_mrr(AllInput{1}, AllInput{2},...
            AllInput{3},AllInput{4},AllInput{5},AllInput{6},AllInput{7},AllInput{8},AllInput{9});
            AllOutput = estr;
        case 14 % mm10_formJ22
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = etau
            [~, ~, ~, ~, ~, ~, ~, ~, ~, etau] = mm10_ehard_mrr(AllInput{1}, AllInput{2},...
            AllInput{3},AllInput{4},AllInput{5},AllInput{6},AllInput{7},AllInput{8},AllInput{9});
            AllOutput = etau;
        case 15 % mm10_formvecs
            % AllInput = props, np1, n, stress, tt, vec1, vec2
            % AllOutput = props, np1, n, stress, tt, vec1, vec2
            [~, ~, ~, ~, ~, vec1, vec2] = mm10_v_mrr(AllInput{1}, AllInput{2},...
            AllInput{3},AllInput{4},AllInput{5},AllInput{6},AllInput{7});
            AllOutput = {vec1, vec2};
        case 16 % mm10_formarrs
            % AllInput = props, np1, n, stress, tt, vec1, vec2, arr1, arr2,both
            % AllOutput = arr1, arr2
        	[~, ~, ~, ~, ~, arr1, arr2] = mm10_a_mrr(AllInput{1}, AllInput{2},...
            AllInput{3}, AllInput{8}, AllInput{9}, AllInput{6}, AllInput{7},AllInput{10});
            AllOutput = {arr1, arr2};
        case 17 % mm10_formvecsi
            error('not a function in matlab')
        case 18 % mm10_output
            % AllInput = props, np1, n, vec1, vec2, np1.stress, np1.tau_tilde
            % AllOutput = np1.slip_incs(1:props.nslip)
            vec1 = AllInput{4};
            slip_incs(1:props.nslip) = vec1(1:props.nslip);
            AllOutput = slip_incs(1:props.nslip);
        case 19 % mm10_tangent->mm10_ed_mrr
            % AllInput = props, np1, n, np1.stress, np1.tau_tilde
            % AllOutput = ed
            [~, ~, ~, ~, ~, ed] = mm10_ed_mrr(AllInput{1}, AllInput{2},...
            AllInput{3}, AllInput{4}, AllInput{5});
            AllOutput = ed;
        case 20 % mm10_tangent->mm10_dgdd_mrr
            % AllInput = props, np1, n, np1.stress, np1.tau_tilde, np1.D
            % AllOutput = dgammadd
            [~, ~, ~, ~, ~, ~, dgammadd] = mm10_dgdd_mrr(AllInput{1}, AllInput{2},...
            AllInput{3}, AllInput{4}, AllInput{5}, AllInput{6});
            AllOutput = dgammadd;
        otherwise
            task
            error('task not implemented')
    end
end

function AllOutput = mm10_plugin_suvic(task,props,AllInput)
% Plugin function for Ma-Roters-Raabe model
    fSchmidt = 0.4082;
%     gamma = (epsdoti/4)/fSchmidt;
%     tau = fSchmidt*sigma;
    switch task
        case 1 % setup_mm10_rknstr
            % AllInput = 0
            % AllOutput = num_hard
            AllOutput = 4;
        case 2 % mm10_init_cc_hist0
            % AllInput = props, history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard),...
            %    history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard)
            % AllOutput = props, history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard), ...
            %    history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard)
            Output1 = AllInput{1};
            Output2 = [props.c1 props.c2 props.c3 props.c4]';
            Output3 = AllInput{3};
            AllOutput = {Output1 Output2 Output3};
        case 3 % mm10_setup
            % AllInput = props, np1, n
            % AllOutput = props, np1, n
            [Output1, Output2, ...
                Output3] = mm10_setup_mrr(AllInput{1}, AllInput{2},AllInput{3});
            AllOutput = {Output1 Output2 Output3};
        case 4 % mm10_formR2
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = h
            np1 = AllInput{2};
            n = AllInput{3};
            stress = AllInput{6};
            s = stress*np1.ms(1:6,8)/fSchmidt; % slip system 8 has the same sign as the applied remote stress
            tt = AllInput{7};
            Bs = tt(1);
            Bl = tt(2);
            R = tt(3);
            K = tt(4);
            locprops = props.p_y; % copy out local properties
            epsdoti = sliprate(s,Bs,Bl,R,K,locprops);
            Bsdot = Bsrate(epsdoti,Bs,locprops);
            Bldot = Blrate(epsdoti,Bl,locprops);
            Rdot = Rrate(epsdoti,R,locprops);
            Kdot = Krate(epsdoti,K,locprops);

            h = zeros(1,props.num_hard);
            dt = np1.tinc;
            h(1) = n.tau_tilde(1) + dt*(Bsdot);
            h(2) = n.tau_tilde(2) + dt*(Bldot);
            h(3) = n.tau_tilde(3) + dt*(Rdot);
            h(4) = n.tau_tilde(4) + dt*(Kdot);
            AllOutput = h;
        case 5 % mm10_formR2i
            error('not a function in matlab')
        case 6 % mm10_form_dbarp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = dbar
%             np1 = AllInput{2};
%             stress = AllInput{6};
%             tt = AllInput{7};
%             dbar = zeros(6,1);
%             Bs = tt(1);
%             Bl = tt(2);
%             R = tt(3);
%             K = tt(4);
%             locprops = props.p_y; % copy out local properties
%             sign8 = sign(stress*np1.ms(1:6,8)/fSchmidt);
%             for i = 1:props.nslip
%                 ms = np1.ms(1:6,i);
%                 rs = stress*ms/fSchmidt; % tau^a
%                 epsdoti = sliprate(abs(rs)*sign8,Bs,Bl,R,K,locprops);
%                 gammadot = (epsdoti/4)/fSchmidt*sign(rs);
%                 dbar = dbar + gammadot*np1.ms(1:6,i);
%             end
%             AllOutput = dbar;
            np1 = AllInput{2};
            vec1 = AllInput{4};
            dbar = np1.ms*vec1(1:props.nslip);
            AllOutput = dbar;
        case 7 % mm10_form_wbarp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = wbar
            np1 = AllInput{2};
            vec1 = AllInput{4};
            wbar = np1.qs*vec1(1:props.nslip);
            AllOutput = wbar;
%         case 7 % mm10_form_wbarp
%             % AllInput = props, np1, n, vec1, vec2, stress, tt
%             % AllOutput = wbar
%             np1 = AllInput{2};
%             stress = AllInput{6};
%             tt = AllInput{7};
%             wbar = zeros(3,1);
%             Bs = tt(1);
%             Bl = tt(2);
%             R = tt(3);
%             K = tt(4);
%             locprops = props.p_y; % copy out local properties
%             sign8 = sign(stress*np1.ms(1:6,8)/fSchmidt);
%             for i = 1:props.nslip
%                 ms = np1.ms(1:6,i);
%                 rs = stress*ms/fSchmidt; % tau^a
%                 epsdoti = sliprate(abs(rs)*sign8,Bs,Bl,R,K,locprops);
%                 gammadot = (epsdoti/4)/fSchmidt*sign(rs);
%                 wbar = wbar + gammadot*np1.qs(1:3,i);
%             end
%             AllOutput = wbar;
        case 8 % mm10_form_wp
            % AllInput = props, np1, n, vec1, vec2, stress, tt
            % AllOutput = w
            np1 = AllInput{2};
            stress = AllInput{6};
            tt = AllInput{7};
            w = zeros(3,1);
            Bs = tt(1);
            Bl = tt(2);
            R = tt(3);
            K = tt(4);
            locprops = props.p_y; % copy out local properties
            sign8 = sign(stress*np1.ms(1:6,8)/fSchmidt);
                rs = stress*np1.ms/fSchmidt; % tau^a
                epsdoti = sliprate(abs(rs)*sign8,Bs,Bl,R,K,locprops);
                gammadot = (epsdoti/4)/fSchmidt.*sign(rs);
                w = np1.qc*gammadot';
%             for i = 1:props.nslip
%                 ms = np1.ms(1:6,i);
%                 rs = stress*ms/fSchmidt; % tau^a
%                 epsdoti = sliprate(abs(rs)*sign8,Bs,Bl,R,K,locprops);
%                 gammadot = (epsdoti/4)/fSchmidt*sign(rs);
%                 w = w + gammadot*np1.qc(1:3,i);
%             end
            AllOutput = w;
        case 9 % mm10_form_dbarpi
            error('not a function in matlab')
        case 10 % mm10_form_wpi
            error('not a function in matlab')
        case 11 % mm10_formJ11
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = dgammadtau
        case 12 % mm10_formJ12
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = dgammadtt
        case 13 % mm10_formJ21
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = estr
        case 14 % mm10_formJ22
            % AllInput = props, np1, n, vec1,vec2,arr1,arr2, stress, tt
            % AllOutput = etau
        case 15 % mm10_formvecs
            % AllInput = props, np1, n, stress, tt, vec1, vec2
            % AllOutput = vec1, vec2
            np1 = AllInput{2};
            stress = AllInput{4};
            tt = AllInput{5};
            Bs = tt(1);
            Bl = tt(2);
            R = tt(3);
            K = tt(4);
            locprops = props.p_y; % copy out local properties
            sign8 = sign(stress*np1.ms(1:6,8)/fSchmidt);
                rs = stress*np1.ms/fSchmidt; % tau^a
                epsdoti = sliprate(abs(rs)*sign8,Bs,Bl,R,K,locprops);
                gammadot = (epsdoti/4)/fSchmidt.*sign(rs);
                AllInput{6} = gammadot';
%             for i = 1:props.nslip
%                 ms = np1.ms(1:6,i);
%                 rs = stress*ms/fSchmidt; % tau^a
%                 epsdoti = sliprate(abs(rs)*sign8,Bs,Bl,R,K,locprops);
%                 gammadot = (epsdoti/4)/fSchmidt*sign(rs);
%                 AllInput{6}(i) = gammadot;
%             end
            AllOutput = {AllInput{6}, AllInput{7}};
        case 16 % mm10_formarrs
            % AllInput = props, np1, n, stress, tt, vec1, vec2, arr1, arr2,both
            % AllOutput = arr1, arr2
        case 17 % mm10_formvecsi
            error('not a function in matlab')
        case 18 % mm10_output
            % AllInput = props, np1, n, vec1, vec2, np1.stress, np1.tau_tilde
            % AllOutput = np1.slip_incs(1:props.nslip)
            vec1 = AllInput{4};
            slip_incs(1:props.nslip) = vec1(1:props.nslip);
            AllOutput = slip_incs(1:props.nslip);
%             slip_incs = zeros(props.nslip,1);
%             np1 = AllInput{2};
%             stress = AllInput{6};
%             tt = AllInput{7};
%             Bs = tt(1);
%             Bl = tt(2);
%             R = tt(3);
%             K = tt(4);
%             locprops = props.p_y; % copy out local properties
%             sign8 = sign(stress*np1.ms(1:6,8)/fSchmidt);
%             for i = 1:props.nslip
%                 ms = np1.ms(1:6,i);
%                 rs = stress*ms/fSchmidt; % tau^a
%                 epsdoti = sliprate(abs(rs)*sign8,Bs,Bl,R,K,locprops);
%                 gammadot = (epsdoti/4)/fSchmidt*sign(rs);
%                 slip_incs(i) = gammadot;
%             end
%             AllOutput = slip_incs;
        case 19 % mm10_tangent->mm10_ed_mrr
            % AllInput = props, np1, n, np1.stress, np1.tau_tilde
            % AllOutput = ed
            AllOutput = zeros(6,props.num_hard);
        case 20 % mm10_tangent->mm10_dgdd_mrr
            % AllInput = props, np1, n, np1.stress, np1.tau_tilde, np1.D
            % AllOutput = dgammadd
            AllOutput = zeros(props.nslip,6);
        otherwise
            task
            error('task not implemented')
    end
end
