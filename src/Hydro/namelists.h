      namelist /CORE/ipre,ibc,ibtp,ntracer_gen,ntracer_age,sed_class,eco_class, &
     &nspool,ihfskip,msc2,mdc2,dt,rnday

      namelist /OPT/ gen_wsett,flag_fib,ics,rearth_pole,rearth_eq,indvel, &
     &imm,ibdef,ihot,ihydraulics,izonal5,slam0,sfea0,iupwind_mom,ihorcon, &
     &hvis_coef0,ishapiro,shapiro0,niter_shap,ihdif,thetai,drampbc, &
     &dramp,nadv,dtb_min,dtb_max,h0,nchi,dzb_min, &
     &hmin_man,ncor,rlatitude,coricoef,nws,wtiminc,iwind_form, &
     &drampwind,iwindoff,ihconsv,isconsv,itur,dfv0,dfh0,h1_pp,h2_pp,vdmax_pp1, &
     &vdmax_pp2,vdmin_pp1,vdmin_pp2,tdmin_pp1,tdmin_pp2,mid,stab,xlsc0, &
     &ibcc_mean,flag_ic,start_year,start_month,start_day,start_hour,utc_start, &
     &itr_met,h_tvd,eps1_tvd_imp,eps2_tvd_imp,ip_weno, &
     &courant_weno,ntd_weno,nquad,epsilon1,epsilon2,epsilon3,ielad_weno,small_elad, &
     &i_prtnftl_weno,inu_tr,step_nu_tr,vnh1,vnh2,vnf1,vnf2, &
     &moitn0,mxitn0,rtol0,iflux,iflux_out_format,inter_mom,h_bcc1,inu_elev,inu_uv, &
     &ihhat,kr_co,rmaxvel,velmin_btrack,btrack_nudge,ibtrack_test,irouse_test, &
     &inunfl,shorewafo,ic_elev,nramp_elev,inv_atm_bnd,prmsl_ref,s1_mxnbt,s2_mxnbt, &
     &iharind,icou_elfe_wwm,drampwafo,nstep_wwm,hmin_radstress,turbinj,turbinjds,alphaw, &
     &fwvor_advxy_stokes,fwvor_advz_stokes,fwvor_gradpress,fwvor_breaking, &
     &fwvor_streaming,fwvor_wveg,fwvor_wveg_NL,wafo_obcramp, &
     &iwbl,cur_wwm,if_source,dramp_ss,ieos_type,ieos_pres,eos_a,eos_b,slr_rate, &
     &rho0,shw,isav,nstep_ice,iunder_deep,h1_bcc,h2_bcc,hw_depth,hw_ratio, &
     &level_age,vclose_surf_frac,iadjust_mass_consv0,ipre2, &
     &ielm_transport,max_subcyc,i_hmin_airsea_ex,hmin_airsea_ex,itransport_only,meth_sink, &
     &iloadtide,loadtide_coef,nu_sum_mult,i_hmin_salt_ex,hmin_salt_ex,h_massconsv,lev_tr_source, &
     &rinflation_icm,iprecip_off_bnd,model_type_pahm
