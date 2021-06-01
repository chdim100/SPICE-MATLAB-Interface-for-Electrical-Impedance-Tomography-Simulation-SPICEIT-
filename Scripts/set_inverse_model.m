function IMDL=set_inverse_model(N,skipcurr,skipvolt,algorithm,filt,lambda,target)
%sets the inverse model where the reconstruction is taking place
%target='d2c2';
%target='d2t3';

imdl2= mk_common_model(target,N);
% Calculate a stimulation pattern
[stim,els] = mk_stim_patterns(N,1,[0,skipcurr+1],[0,skipvolt+1],{},1);
% Solve all voltage patterns
imdl2.fwd_model.stimulation = stim;
imdl2.fwd_model.meas_select=els;
imdl2.fwd_solve.get_all_meas = 1;
%set algorithm parameters
imdl2.hyperparameter.value = lambda;
switch algorithm
    case{1,'direct','one_step'}
        imdl2.solve=@inv_solve_diff_GN_one_step;
    case{2,'iterational','mulistep','iteration','GN'}
        imdl2.solve=@inv_solve_abs_GN;
        imdl2.inv_solve_gn.max_iterations = 20;
    case{3,'TV','Total_Variation','Pdipm','pdipm','PDIPM'}
        imdl2.solve=@inv_solve_TV_pdipm;
        imdl2.R_prior=     @prior_TV;
        imdl2.parameters.term_tolerance= 1e-3;
        imdl2.parameters.max_iterations=25;
        imdl2.parameters.keep_iterations= 25;
    case{4,'Back-Projection','Back-projection','BP'}
        imdl2.solve= @inv_solve_backproj;
        imdl2.inv_solve_backproj.type='simple_filter';
    case{5,'iterationalGN','mulistepGN','iterationGN','GN2'}
        imdl2.solve=@inv_solve_gn;
        imdl2.inv_solve_gn.max_iterations = 20;
        imdl2.inv_solve_gn.return_working_variables = 1;
    otherwise
        error('under construction')
end
if algorithm~=4
    switch filt
        case {1,'NOSER','Noser','noser'}
            imdl2.RtR_prior= @prior_noser;
        case{2,'LAPLACE','Laplace','laplace'}
            imdl2.RtR_prior= @prior_laplace;
        case{3,'Gaussian','Gauss','GAUSSIAN','gaussian'}
            imdl2.RtR_prior = @prior_gaussian_HPF;
        case{4,'TV','Total_Variation'}
            imdl2.RtR_prior=@prior_TV;
        otherwise
            error('under construction')
    end
end
imdl2.jacobian_bkgnd.value= 1;
IMDL=imdl2;
end