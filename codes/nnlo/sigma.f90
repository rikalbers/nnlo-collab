! Source containing routines related to the organization of the 
! cross section calculation:
subroutine NNLOcalc()
use flags
use statistics
use histo
use phasespace
use scales
implicit none
!
!
!
! Initialization:
  call init_stat
  if (flg_analysis) call init_analysis
!
! Integration:
! The integration is two-fold hidden in order to 
! be able to introduce new integrators in a relatively easy way:
  call integrate
!
! We finalize the statistics and write out the histos:
  call output_stat
  if (flg_analysis) call output_hist
!
end subroutine NNLOcalc
!
subroutine CalcSigma(cont)
use flags
use subprocesses
use phasespace
use statistics
use histo
use my_scale
use coupling
use scales
use regions
implicit none
!
  character (len=5) , intent(in) :: cont
!
  integer :: iproc,icut
  logical :: discard
  real(kind(1d0)) :: cfunc
  real(kind(1d0)) :: weight,weightPS
  real(kind(1d0)) :: weight_R,weight_R_A1
!
  real(kind(1d0)) :: minsij,smeRV,smeR,subCir,subSr,subCirSr
!
  interface
    subroutine sigmaB(p,weightPS,dsigmaBi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaBi
    !
    end subroutine sigmaB
!
    subroutine sigmaV(p,weightPS,dsigmaVi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaVi
    !
    end subroutine sigmaV
!
    subroutine sigmaR(p,weightPS,dsigmaRi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRi
    !
    end subroutine sigmaR
!
    subroutine sigmaRV(p,weightPS,dsigmaRVi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRVi
    !
    end subroutine sigmaRV
!
    subroutine sigmaRR(p,weightPS,dsigmaRRi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRRi
    !
    end subroutine sigmaRR
!
    subroutine sigmaR_A1(p,weightPS,Cir,Sr,CSir,dsigmaR_A1)
    use particles
    use regions
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CSir
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaR_A1
    !
    end subroutine sigmaR_A1
!
    subroutine sigmaRR_A1(p,weightPS,Cir,Sr,CSir,dsigmaRR_A1)
    use particles
    use regions
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CSir
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRR_A1
    !
    end subroutine sigmaRR_A1
!
    subroutine sigmaR_I1(p,weightPS,dsigmaI1)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaI1
!
    end subroutine sigmaR_I1
!
    subroutine sigma_B_Is(p,Bij,weightPS,dsigmaB,dsigmaV,dsigmaVV)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaB,  &
                                                      dsigmaV,  &
                                                      dsigmaVV
!
    end subroutine sigma_B_Is
!
    subroutine sigma_B_V_Is(p,Bij,BijLaurent,Bijkl,BijklLaurent, &
                            Vij,VijLaurent,weightPS, &
                            dsigmaB,dsigmaV,dsigmaVV)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij,Vij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
      real(kind(1d0)) , dimension(:,:,:,:,-4:) , intent(out) :: BijklLaurent
      real(kind(1d0)) , dimension(:,:,-4:) , intent(out) :: BijLaurent
      real(kind(1d0)) , dimension(:,:,-4:) , intent(out) :: VijLaurent
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaB,  &
                                                      dsigmaV,  &
                                                      dsigmaVV
!
    end subroutine sigma_B_V_Is
!
    subroutine sigma_R_Is(p,Rij,weightPS,dsigmaR,dsigmaRV)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaR, &
                                                      dsigmaRV
!
    end subroutine sigma_R_Is
!
    subroutine Calc_mp1_subs(p,ptilde,Bij,Vij,Bijk,Bijkl,Bmunuij, &
                             weightPS, &
                             Cir,Sr,CirSr, &
                             subsR,subsRV,subsRRA1)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , dimension(:,:,:) , intent(out) :: Bijk
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
      real(kind(1d0)) , intent(inout) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CirSr
      real(kind(1d0)) , dimension(:) , intent(out) :: subsR
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRV
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRRA1
!
    end subroutine Calc_mp1_subs
!
    subroutine Calc_mp2_subs(p,phat,ptilde,         &
                             Bij,Bijkl,Bmunuij,Rij, &
                             weightPS,              &
                             Cir,Sr,CirSr,          &
                             Cirs,CSirs,Cirjs,Srs,  &
                             subsRR)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CirSr
      type(subterms) , intent(in) :: Cirs,CSirs,Cirjs,Srs
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRR
!
    end subroutine Calc_mp2_subs
!
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine analysis(p,weights)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:) , intent(in) :: weights
    !
    end subroutine analysis
  end interface
!
! We can have several different multiplicities
! Born kinematics:
  if (cont.eq.'Bkin ') then
! Nullifying arrays holding the Born and regularized virtual
! contributions, arrays are defined due to scale uncertainty 
! studies:
    dsigmaB  = 0d0
    dsigmaV  = 0d0
    dsigmaVV = 0d0
    do iproc=1,num_flv_ch_Bkin
      if (flg_LO) dsigmaBi = 0d0
      if (flg_NLO_V) dsigmaVi = 0d0
      if (flg_NNLO_VV) dsigmaVVi = 0d0
! Obtain the phase space point for each subprocess:
      call GetPSpoint(iproc,cont,discard,weightPS)
! Associate flavor information to momenta:
      call CreateParts(ps_B,flv_ch_Bkin(:,iproc),parts_B)
! Multiply with the weight coming from subprocesses a-like:
      weightPS = weightPS *  weights_Bkin(iproc)
! Apply cuts on the PS point:
      call Cuts(parts_B,cfunc)
      if (cfunc.ne.0d0) then
! Old way of doing the VV line:
! Calculation of the Born part and all the I operators:
!        call sigma_B_Is(parts_B,Bij_arr,weightPS, &
!                        dsigmaBi,dsigmaB_I1,dsigmaB_I2)
!        call sigmaV(parts_B,weightPS,dsigmaVi)
!        if (flg_NLO_V) then
!          dsigmaV = dsigmaV + (dsigmaVi + dsigmaB_I1) * cfunc
!        end if
!
! Obtaining the Born, the virtual and all the integrated 
! counterterms in one routine to optimize the call to the SMEs:
! dsigmaBi only contains the Born weight, while dsigmaVi contains
! not just the virtual but the appropriate NLO integrated 
! counterterms too,
! finally dsigma_Innlo contains all the NNLO-related integrated
! counterterms.
        call sigma_B_V_Is(parts_B, &
                          Bij_arr,BijLaurent_arr, &
                          Bijkl_arr,BijklLaurent_arr, &
                          Vij_arr,VijLaurent_arr, &
                          weightPS, &
                          dsigmaBi,dsigmaVi,dsigmaVVi)
        if (flg_LO) dsigmaB = dsigmaB + dsigmaBi * cfunc
! dsigmaVi already contains the the integrated counterterms needed for
! the pure NLO calculation:
        if (flg_NLO_V) then
          dsigmaV = dsigmaV + dsigmaVi * cfunc
        end if
! For the NNLO VV part we have to add the VV contribution and
! the corresponding integrated subtraction terms:
        if (flg_NNLO_VV) then
          dsigmaVV = dsigmaVV + dsigmaVVi * cfunc
        end if
!
      end if
!
      weight = (dsigmaBi(centscale)   &
             +  dsigmaVi(centscale)   &
             +  dsigmaVVi(centscale)  &
                ) * cfunc
! This routine gives back the weight to the integrator:
      call PutWeight(iproc,cont,weight)
!
! If needed the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(parts_B,dsigmaBi + dsigmaVi + dsigmaVVi)
      end if
!
    end do
!
    if (flg_LO) call accu_stat('B ',dsigmaB)
    if (flg_NLO_V) call accu_stat('V ',dsigmaV)
    if (flg_NNLO_VV) call accu_stat('VV',dsigmaVV)
    if (.not.flg_optim.and.flg_analysis) call accu_hist
!
  elseif (cont.eq.'Rkin ') then
    dsigmaR = 0d0
    dsigmaRV = 0d0
    do iproc=1,num_flv_ch_Rkin
      if (flg_NLO_R) dsigmaRi = 0d0
      if (flg_NNLO_RV) dsigmaRVi = 0d0
      if (flg_NNLO_RV) dsigmaR_I1 = 0d0
! We obtain the phase space point for each subprocess:
      call GetPSpoint(iproc,cont,discard,weightPS)
! We associate flavor information to momenta:
      call CreateParts(ps_R,flv_ch_Rkin(:,iproc),parts_R)
! We have to deal with all the others too:
      weightPS = weightPS * weights_Rkin(iproc)
! Start with calculating the subtractions:
! The reason behind this ordering is that it is possible that
! the Real phase space passed the ymin cut but the underlying
! Born one does not resulting in a veto on the PS point, hence
! it is useless to calculate the Real and/or the Real-Virtual:
!
! Master routine manages all the subtractions defined for
! m+1 partonic final states:
      call Calc_mp1_subs(parts_R,parts_B,                        &
                         Bij_arr,Vij_arr,                        &
                         Bijk_arr,Bijkl_arr,Bmunuij_arr,         &
                         weightPS,                               &
                         subterms_Cir_R(iproc),                  &
                         subterms_Sr_R(iproc),                   &
                         subterms_CSir_R(iproc),                 &
                         subs_R,subs_RV,subs_RRA1)
! R and RV is only calculated with the weight is nonzero:
      if (weightPS.ne.0) then
! Apply cuts on the PS point:
        call Cuts(parts_R,cfunc)
        if (cfunc.ne.0d0) then
! Calculation of the Real part and the I operator:
          call sigma_R_Is(parts_R,Rij_arr,weightPS, &
                          dsigmaRi,dsigmaR_I1)
          if (flg_NLO_R) dsigmaR = dsigmaR + dsigmaRi * cfunc
!
          call sigmaRV(parts_R,weightPS,dsigmaRVi)
          if (flg_NNLO_RV) dsigmaRV = dsigmaRV &
                                    + (dsigmaRVi + dsigmaR_I1) * cfunc
        end if
!
! Weight for the PS generator is constructed from R, RV and subtractions
! R, RV and I_1^{(0)}\times R receives an additional weight coming from the
! cut function:
        weight = (dsigmaRi(centscale)    &
               +  dsigmaRVi(centscale)   &
               +  dsigmaR_I1(centscale)) &
               * cfunc                   &
               -  subs_R(centscale)      &
               -  subs_RV(centscale)     &
               -  subs_RRA1(centscale)
!
! Constructing entries for the statistics:
        if (flg_NLO_R.and.flg_NLO_R_A1) dsigmaR = dsigmaR - subs_R
        if (flg_NNLO_RV.and.flg_NNLO_RV_A1) dsigmaRV = dsigmaRV - subs_RV
        if (flg_NNLO_RV.and.flg_NNLO_RRA1_A1) dsigmaRV = dsigmaRV - subs_RRA1
!
! If needed the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then 
          call analysis(parts_R,dsigmaRi + dsigmaRVi + dsigmaR_I1)
        end if
! If the weight of the PS point became zero due to an underlying Born
! PS point not passing the ymin cut the weight for the integrator should be
! zero:
      else
        weight = 0
      end if
! We have to give back a weight to the PS generator/integrator,
! we always choose the central scale choice:
      call PutWeight(iproc,cont,weight)
!
    end do
    if (flg_NLO_R) call accu_stat('R ',dsigmaR)
    if (flg_NNLO_RV) call accu_stat('RV',dsigmaRV)
    if (.not.flg_optim.and.flg_analysis) call accu_hist
!
  elseif (cont.eq.'RRkin') then
    dsigmaRR = 0d0
    do iproc=1,num_flv_ch_RRkin
      if (flg_NNLO_RR) dsigmaRRi = 0d0
! We obtain the phase space point for each subprocess:
      call GetPSpoint(iproc,cont,discard,weightPS)
! We associate flavor information to momenta:
      call CreateParts(ps_RR,flv_ch_RRkin(:,iproc),parts_RR)
! We have to deal with all the others too:
      weightPS = weightPS * weights_RRkin(iproc)
! Apply cuts on the PS point:
      call Cuts(parts_RR,cfunc)
      if (cfunc.ne.0d0) then
        call sigmaRR(parts_RR,weightPS,dsigmaRRi)
        if (flg_NNLO_RR) dsigmaRR = dsigmaRR + dsigmaRRi * cfunc
!
      end if
!
! Calling the master routine to calculate all subtraction terms
! defined for the RR term m+2 partonic final state:
      call Calc_mp2_subs(parts_RR,parts_R,parts_B, &
                         Bij_arr,Bijkl_arr,        &
                         Bmunuij_arr,Rij_arr,      &
                         weightPS,                 &
                         subterms_Cir_RR(iproc),   &
                         subterms_Sr_RR(iproc),    &
                         subterms_CSir_RR(iproc),  &
                         subterms_Cirs_RR(iproc),  &
                         subterms_CSirs_RR(iproc), &
                         subterms_Cirjs_RR(iproc), &
                         subterms_Srs_RR(iproc),   &
                         subs_RR)
!
      weight = dsigmaRRi(centscale) * cfunc &
             -  subs_RR(centscale)
! We have to give back a weight to the PS generator/integrator,
! we always choose the central scale choice:
      call PutWeight(iproc,cont,weight)
      if (flg_NNLO_RR.and. &
          (flg_NNLO_RR_A1.or.flg_NNLO_RR_A2.or.flg_NNLO_RR_A12)) then
        dsigmaRR = dsigmaRR - subs_RR
      end if
!
! If needed the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(parts_RR,dsigmaRRi)
      end if
!
!      print *,"weight: ",weight
!
    end do
    if (flg_NNLO_RR) call accu_stat('RR',dsigmaRR)
    if (.not.flg_optim.and.flg_analysis) call accu_hist
  end if
!
end subroutine CalcSigma
!
subroutine sigmaB(p,weightPS,dsigmaBi)
use particles
use flags
use subprocesses
use statistics
use my_scale
use coupling
use scales
use misc
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaBi
!
  integer :: iscale
  real(kind(1d0)) :: sme
  real(kind(1d0)) :: fluxfact
!
  interface
    subroutine CalcB(parts,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  dsigmaBi = 0d0
!
! If the weight is zero or contribution is not needed we do nothing:
  if ((weightPS.eq.0d0).or. &
      (.not.flg_LO)) then
    return
  end if
!
! We can have a dynamical scale hence we make a call to the
! scale routine:
  call calcmyscales(p)
! For the Born part no scale dependence can come from the SME,
! hence it is only called once:
  call CalcB(p,sme)
! We include the couplings:
  call calc_couplings
! We calculate the flux factor: 
  fluxfact = calcFluxFact(p)
! We store the original xiR and xiF values:
  call StoreXis_scales
! Evaluating possible scales:
  do iscale=1,nscales
! We change the scales and recalculate the couplings accordingly:
    call ChangeXis_scales(iscale)
! We recalculate the couplings:
    call calcmyscales(p)
    call calc_couplings
! We store the weight:
!        print *,"xir,xif: ",xir,xif
    dsigmaBi(iscale) = weightPS * sme * fluxfact &
                     * alphas**border_as * alphaEM**border_aEM
! Next line is only for testing with g_s = e = 1:
!                     / (4d0*pi)**(border_as+border_aEM)
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
  return
!
end subroutine sigmaB
!
subroutine sigmaV(p,weightPS,dsigmaVi)
use particles
use flags
use subprocesses
use statistics
use my_scale
use coupling
use scales
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaVi
!
  integer :: iscale
  real(kind(1d0)) :: sme
  real(kind(1d0)) :: fluxfact
!
  interface
    subroutine CalcV(parts,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeVLaurent
!
    end subroutine CalcV
  end interface
!
  dsigmaVi = 0d0
!
! If the weight is zero we do not anything:
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NLO_V)) then
    return
  end if
!
! We can have a dynamical scale hence we make a call to the
! scale routine:
  call calcmyscales(p)
! We calculate the flux factor: 
  fluxfact = calcFluxFact(p)
! We store the original xiR and xiF values:
  call StoreXis_scales
! Evaluating possible scales:
  do iscale=1,nscales
! We change the scales and recalculate the couplings accordingly:
    call ChangeXis_scales(iscale)
! We recalculate the couplings:
    call calcmyscales(p)
    call calc_couplings
    call CalcV(p,sme)
! We store the weight:
!    print *,"xir,xif: ",xir,xif
    dsigmaVi(iscale) = sme &
                     * weightPS &
                     * fluxfact &
                     * alphas**(border_as+1) &
                     * alphaEM**border_aEM
! Next line is only for testing with g_s = e = 1:
!                     / (4d0*pi)**(border_as+1+border_aEM)
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
  return
!
end subroutine sigmaV
!
subroutine sigmaR(p,weightPS,dsigmaRi)
use particles
use flags
use statistics
use my_scale
use coupling
use scales
use misc
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRi
!
  integer :: iscale
  real(kind(1d0)) :: sme
  real(kind(1d0)) :: fluxfact
!
  interface
    subroutine CalcR(parts,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
  end interface
!
  dsigmaRi = 0d0
!
! If the weight is zero we do not anything:
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NLO_R)) then
    return
  end if
!
! We can have a dynamical scale hence we make a call to the
! scale routine:
  call calcmyscales(p)
! For the Real part no scale dependence can come from the SME,
! if the couplings are factored out hence it is only called once:
  call CalcR(p,sme)
! We include the couplings:
  call calc_couplings
! We calculate the flux factor:
  fluxfact = CalcFluxFact(p)
! We store the original xiR and xiF values:
  call StoreXis_scales
! Evaluating possible scales:
  do iscale=1,nscales
! We change the scales and recalculate the couplings accordingly:
    call ChangeXis_scales(iscale)
! We recalculate the couplings:
    call calcmyscales(p)
    call calc_couplings
! We store the weight:
!    print *,"xir,xif: ",xir,xif
    dsigmaRi(iscale) = weightPS * sme * fluxfact &
                     * alphas**(border_as+1) * alphaEM**border_aEM
! Next line is only for testing with g_s = e = 1:
!                     / (4d0*pi)**(border_as+1+border_aEM)
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
  return
!
end subroutine sigmaR
!
subroutine sigmaRR(p,weightPS,dsigmaRRi)
use particles
use flags
use statistics
use my_scale
use coupling
use scales
use misc
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRRi
!
  integer :: icut
  integer :: iscale
  real(kind(1d0)) :: sme,cfunc,weight
  real(kind(1d0)) :: fluxfact
!
  interface
    subroutine apply_cuts(p,icut)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      integer , intent(out) :: icut
!
    end subroutine apply_cuts
!
    subroutine cutfunc(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine cutfunc
!
    subroutine analysis(p,weights)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:) , intent(in) :: weights
    !
    end subroutine analysis
!
    subroutine CalcRR(parts,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
  end interface
!
  dsigmaRRi = 0d0
!
! If the weight is zero we do not anything:
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NNLO_RR)) then
    return
  end if
!
  weight = weightPS
!
! We can have a dynamical scale hence we make a call to the
! scale routine:
  call calcmyscales(p)
! For the Real part no scale dependence can come from the SME,
! if the couplings are factored out hence it is only called once:
  call CalcRR(p,sme)
! We include the couplings:
  call calc_couplings
! We calculate the flux factor:
  fluxfact = CalcFluxFact(p)
! We store the original xiR and xiF values:
  call StoreXis_scales
! Evaluating possible scales:
  do iscale=1,nscales
! We change the scales and recalculate the couplings accordingly:
    call ChangeXis_scales(iscale)
! We recalculate the couplings:
    call calcmyscales(p)
    call calc_couplings
! We store the weight:
!    print *,"xir,xif: ",xir,xif
    dsigmaRRi(iscale) = weight * fluxfact * sme &
                      * alphas**(border_as+2) * alphaEM**border_aEM
! Next line is only for testing with g_s = e = 1:
!                      / (4d0*pi)**(border_as+2+border_aEM)
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
  return
!
end subroutine sigmaRR
!
subroutine sigmaR_A1(p,weightPS,Cir,Sr,CSir,dsigmaA1)
use regions
use particles
use flags
use phasespace
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  type(subterms) , intent(in) :: Cir,Sr,CSir
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaA1
!
!
!
  interface
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcA1subtraction(p,ptilde,Bij,weightPS,Cir,Sr,CSir, & 
                                 CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
                                 A1term)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CSir
      real(kind(1d0)) , dimension(:) , intent(out) :: A1term
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine CalcA1subtraction
!
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
  end interface
!
  dsigmaR_A1 = 0d0
!
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NLO_R_A1)) then
    return
  end if
! We calculate the invariants:
  call CalcSubInvariants(parts_R)
! Calculation of the subtraction terms:
  call CalcA1subtraction(parts_R,parts_B,Bij_arr,weightPS, &
                         Cir,Sr,CSir, &
                         CalcB,CalcBmunu,CalcBij,dsigmaA1)
!
end subroutine sigmaR_A1
!
subroutine sigmaB_I1(p,weightPS,dsigmaI1)
use particles
use flags
use phasespace
use my_scale
use coupling
use scales
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaI1
!
  integer :: iscale
  real(kind(1d0)) :: fluxfact,I1term
  real(kind(1d0)) :: smeB
!
!
  interface
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcI1nlo(p,smeB,Bij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: smeB
      real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
    end subroutine CalcI1nlo
  end interface
!
  dsigmaI1 = 0d0
!
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NLO_V)) then
    return
  end if
!
! We can have a dynamical scale hence we make a call to the
! scale routine:
  call calcmyscales(p)
! Include the couplings:
  call calc_couplings
! Calculate the flux factor: 
  fluxfact = calcFluxFact(p)
! Calculating the color-correlated Born SME:
  call CalcBij(p,Bij_arr)
! Obtaining the SME from the above result:
  call CastSMEijToSME(p,Bij_arr,smeB)
! Store the original xiR and xiF values:
  call StoreXis_scales
! Evaluating possible scales:
  do iscale=1,nscales
! Change the scales and recalculate the couplings accordingly:
    call ChangeXis_scales(iscale)
! Recalculate the couplings:
    call calcmyscales(p)
    call calc_couplings
    call CalcI1nlo(p,smeB,Bij_arr,I1term)
!
    dsigmaI1(iscale) = weightPS * I1term * fluxfact &
                     * alphas**(border_as+1) * alphaEM**border_aEM
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
  return
!
end subroutine sigmaB_I1
!
subroutine Cuts(p,cfunc)
use flags
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(out) :: cfunc
!
  integer :: icut
!
!
  interface
    subroutine apply_cuts(p,icut)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      integer , intent(out) :: icut
!
    end subroutine apply_cuts
!
    subroutine cutfunc(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine cutfunc
  end interface
!
! By default the output is zero which means not passing the
! the cuts:
  cfunc = 0d0
! If cuts are requested, call the cuts routine:
  if (flg_cuts) then
    call apply_cuts(p,icut)
    cfunc = cfunc + icut
  end if
! Can have a cut function, too:
  if (flg_cutfunc) call cutfunc(p,cfunc)
!
! If no cuts are present giving back one:
  if ((.not.flg_cuts).and.(.not.flg_cutfunc)) cfunc = 1d0
!
  return
!
end subroutine Cuts
!
subroutine sigmaRV(p,weightPS,dsigmaRVi)
use particles
use flags
use subprocesses
use statistics
use my_scale
use coupling
use scales
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRVi
!
  integer :: iscale
  real(kind(1d0)) :: sme
  real(kind(1d0)) :: fluxfact
!
  interface
    subroutine CalcRV(parts,smeRV,smeRVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeRVLaurent
!
    end subroutine CalcRV
  end interface
!
  dsigmaRVi = 0d0
!
! If the weight is zero we do not anything:
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NNLO_RV)) then
    return
  end if
!
! We can have a dynamical scale hence we make a call to the
! scale routine:
  call calcmyscales(p)
! We calculate the flux factor: 
  fluxfact = calcFluxFact(p)
! We store the original xiR and xiF values:
  call StoreXis_scales
! Evaluating possible scales:
  do iscale=1,nscales
! We change the scales and recalculate the couplings accordingly:
    call ChangeXis_scales(iscale)
! We recalculate the couplings:
    call calcmyscales(p)
    call calc_couplings
    call CalcRV(p,sme)
! We store the weight:
!    print *,"xir,xif: ",xir,xif
    dsigmaRVi(iscale) = weightPS * sme * fluxfact &
                      * alphas**(border_as+2) * alphaEM**border_aEM
! Next line is only for testing with g_s = e = 1:
!                      / (4d0*pi)**(border_as+2+border_aEM)
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
  return
!
end subroutine sigmaRV
!
subroutine sigmaR_I1(p,weightPS,dsigmaI1)
use particles
use flags
use phasespace
use my_scale
use coupling
use scales
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaI1
!
  integer :: iscale
  real(kind(1d0)) :: fluxfact,I1term
! TEST:
  real(kind(1d0)) :: Virt
  real(kind(1d0)) , dimension(-4:2) :: I1Laurent,VirtLaurent
!
!
  interface
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine CalcI1(p,Bij,CalcSMEB,CalcSMEBij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
    end subroutine CalcI1
  end interface
!
  dsigmaI1 = 0d0
!
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NNLO_R_I1)) then
    return
  end if
!
! We can have a dynamical scale hence we make a call to the
! scale routine:
  call calcmyscales(p)
! Include the couplings:
  call calc_couplings
! Calculate the flux factor: 
  fluxfact = calcFluxFact(p)
! Store the original xiR and xiF values:
  call StoreXis_scales
! Evaluating possible scales:
  do iscale=1,nscales
! Change the scales and recalculate the couplings accordingly:
    call ChangeXis_scales(iscale)
! Recalculate the couplings:
    call calcmyscales(p)
    call calc_couplings
    call CalcI1(p,Rij_arr,CalcR,CalcRij,I1term)
!
    dsigmaI1(iscale) = weightPS * I1term * fluxfact &
                     * alphas**(border_as+2) * alphaEM**border_aEM
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
  return
!
end subroutine sigmaR_I1
!
subroutine sigmaRR_A1(p,weightPS,Cir,Sr,CSir,dsigmaA1)
use regions
use particles
use flags
use phasespace
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: weightPS
  type(subterms) , intent(in) :: Cir,Sr,CSir
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaA1
!
!
!
  interface
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcA1subtraction(p,ptilde,Bij,weightPS,Cir,Sr,CSir, & 
                                 CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
                                 A1term)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CSir
      real(kind(1d0)) , dimension(:) , intent(out) :: A1term
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine CalcA1subtraction
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine CalcRmunu
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
  end interface
!
  dsigmaRR_A1 = 0d0
!
  if ((weightPS.eq.0d0).or. &
      (.not.flg_NNLO_RR_A1)) then
    return
  end if
! We calculate the invariants:
  call CalcSubInvariants(parts_RR)
! Calculation of the subtraction terms:
  call CalcA1subtraction(parts_RR,parts_R,Rij_arr,weightPS, &
                         Cir,Sr,CSir, &
                         CalcR,CalcRmunu,CalcRij,dsigmaA1)
!
end subroutine sigmaRR_A1
!
! This routine calculates the Born and all the I operators:
subroutine sigma_B_Is(p,Bij,weightPS,dsigmaB,dsigmaV,dsigmaVV)
use particles
use flags
use scales
use my_scale
use coupling
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaB,  &
                                                  dsigmaV,  &
                                                  dsigmaVV
!
  integer :: iscale
  real(kind(1d0)) :: smeB,I1
  real(kind(1d0)) :: fluxfact
!
!
  interface
    subroutine CalcB(parts,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcI1nlo(p,smeB,Bij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: smeB
      real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
    end subroutine CalcI1nlo
  end interface
!
  dsigmaB  = 0d0
  dsigmaV  = 0d0
  dsigmaVV = 0d0
!
! if the weight is zero or both of the contributions are
! not needed return back:
  if ((weightPS.eq.0).or. &
      ((.not.flg_LO).and. &
       (.not.flg_NLO_V).and. &
       (.not.flg_NNLO_VV))) then
    return
  end if
!
! Calculation of the flux factor:
  fluxfact = calcFluxFact(p)
! Dealing only with tree-level amplitudes hence no scale
! dependence occurs if couplings are factored out.
! If only the Born is needed the Born SME is calculated only:
  if (flg_LO.and..not.(flg_NLO_V.or.flg_NNLO_VV)) then
    call calcB(p,smeB)
! If the virtual is needed the corresponding I operator needs
! the color-correlated Born SMEs, hence the Born SME is obtained from
! the color-correlated one:
  elseif (flg_NLO_V.and..not.flg_NNLO_VV) then
    call calcBij(p,Bij)
    call CastSMEijToSME(p,Bij,smeB)
  elseif (flg_NNLO_VV) then
    print *,"The VV part is not yet implemented"
    stop "sigmaB_B_Is"
  end if
! Storing original xi values:
  call StoreXis_scales
! Run through all the scales:
  do iscale=1,nscales
! Changing the scales and recalculating the couplings accordingly:
    call ChangeXis_scales(iscale)
! Setting up the new muR and muF:
    call calcmyscales(p)
! Recalculating the couplings:
    call calc_couplings
    if (flg_LO) then
      dsigmaB(iscale) = smeB &
                      * weightPS &
                      * fluxfact &
                      * alphas**border_as &
                      * alphaEM**border_aEM
    end if
    if (flg_NLO_V) then
! Obtaining the I1 operator:
      call CalcI1nlo(p,smeB,Bij,I1)
      dsigmaV(iscale) = I1 &
                      * weightPS &
                      * fluxfact &
                      * alphas**(border_as+1) &
                      * alphaEM**border_aEM
    end if
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
!
end subroutine sigma_B_Is
!
! This routine calculates the Real and the corresponding I operator
! coming from the RR_A1:
subroutine sigma_R_Is(p,Rij,weightPS,dsigmaR,dsigmaRV)
use particles
use flags
use scales
use my_scale
use coupling
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaR, &
                                                  dsigmaRV
!
  integer :: iscale
  real(kind(1d0)) :: smeR,I1
  real(kind(1d0)) :: fluxfact
!
!
  interface
    subroutine CalcR(parts,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine CalcI1nlo(p,smeB,Bij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: smeB
      real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
    end subroutine CalcI1nlo
  end interface
!
  dsigmaR  = 0d0
  dsigmaRV = 0d0
!
! if the weight is zero or both of the contributions are
! not needed return back:
  if ((weightPS.eq.0).or. &
      ((.not.flg_NLO_R).and. &
       (.not.flg_NNLO_RV))) then
    return
  end if
!
! Calculation of the flux factor:
  fluxfact = calcFluxFact(p)
! Dealing only with tree-level amplitudes hence no scale
! dependence occurs if couplings are factored out.
! If only the Real is needed the Real SME is calculated only:
  if (flg_NLO_R.and..not.flg_NNLO_RV) then
    call calcR(p,smeR)
! If the real-virtual is needed the corresponding I operator needs
! the color-correlated Real SMEs, hence the Real SME is obtained from
! the color-correlated one:
  elseif (flg_NNLO_RV) then
    call calcRij(p,Rij)
    call CastSMEijToSME(p,Rij,smeR)
  end if
! Storing original xi values:
  call StoreXis_scales
! Run through all the scales:
  do iscale=1,nscales
! Changing the scales and recalculating the couplings accordingly:
    call ChangeXis_scales(iscale)
! Setting up the new muR and muF:
    call calcmyscales(p)
! Recalculating the couplings:
    call calc_couplings
    if (flg_NLO_R) then
      dsigmaR(iscale) = smeR &
                      * weightPS &
                      * fluxfact &
                      * alphas**(border_as+1) &
                      * alphaEM**border_aEM
    end if
    if (flg_NNLO_RV) then
! Obtaining the I1 operator:
      call CalcI1nlo(p,smeR,Rij,I1)
      dsigmaRV(iscale) = I1 &
                       * weightPS &
                       * fluxfact &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
    end if
  end do  
! Restoring the old xi values:
  call RestoreXis_scales
!
!
end subroutine sigma_R_Is
!
! Routine to calculate the Born, the virtual ( + the NLO 
! integrated counterterms) and the NNLO integrated counterterms:
subroutine sigma_B_V_Is(p,Bij,BijLaurent,Bijkl,BijklLaurent, &
                        Vij,VijLaurent,weightPS, &
                        dsigmaB,dsigmaV,dsigmaVV)
use particles
use flags
use scales
use my_scale
use coupling
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij,Vij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:,-4:) , intent(out) :: BijklLaurent
  real(kind(1d0)) , dimension(:,:,-4:) , intent(out) :: BijLaurent
  real(kind(1d0)) , dimension(:,:,-4:) , intent(out) :: VijLaurent
  real(kind(1d0)) , intent(in) :: weightPS
  real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaB,  &
                                                  dsigmaV,  &
                                                  dsigmaVV
!
  integer :: iscale
  real(kind(1d0)) :: smeB,smeV,smeVV,I_NLO,I_NNLO
  real(kind(1d0)) :: fluxfact
  real(kind(1d0)) , dimension(-4:2) :: I_NLOLaurent
!
!DEBUG!!!!
  integer :: i
  real(kind(1d0)) , dimension(-4:2) :: BornLaurent
  real(kind(1d0)) , dimension(-4:2) :: Laurent1,Laurent2
  real(kind(1d0)) :: I_NLO1,I_NLO2,V1,V2,temp,B1,B2
!
!
  interface
    subroutine CalcB(parts,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
!
    subroutine CalcV(parts,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeVLaurent
!
    end subroutine CalcV
!
    subroutine CalcVddim(parts,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeVLaurent
!
    end subroutine CalcVddim
!
    subroutine CalcVij(p,Vij,VijLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine CalcVij
!
    subroutine CalcVijddim_new(p,Vij,VijLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine CalcVijddim_new
!
    subroutine CalcVV(parts,smeVV,smeVVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeVV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeVVLaurent
!
    end subroutine CalcVV
!
    subroutine CalcI1nlo(p,smeB,Bij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: smeB
      real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
    end subroutine CalcI1nlo
!
    subroutine CalcInnlo(p,                  &
                         smeB,BornLaurent,   &
                         Bij,BijLaurent,     &
                         Bijkl,BijklLaurent, &
                         smeV,VijLaurent,    &
                         Inloterm,Innloterm,InloLaurent,InnloLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: smeB,smeV
      real(kind(1d0)) , dimension(-4:2) :: BornLaurent
      real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
      real(kind(1d0)) , dimension(:,:,-4:) :: BijLaurent
      real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: Bijkl
      real(kind(1d0)) , dimension(:,:,:,:,-4:) :: BijklLaurent
      real(kind(1d0)) , dimension(:,:,-4:) , intent(in) :: VijLaurent
      real(kind(1d0)) , intent(out) :: Inloterm,Innloterm
      real(kind(1d0)) , dimension(-4:2) , intent(out) :: InloLaurent
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: InnloLaurent
!
    end subroutine CalcInnlo
!
! DEBUG!!!!
    subroutine CalcBijddim(p,Bij,BijLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: &
        BijLaurent
!
    end subroutine CalcBijddim
!
    subroutine CalcBijklddim(p,Bijkl,BijklLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
      real(kind(1d0)) , optional , dimension(:,:,:,:,-4:) , intent(out) :: &
        BijklLaurent
!
    end subroutine CalcBijklddim
!
    subroutine CalcI1nloddim(p,BLaurent,BijLaurent,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(-4:2) , intent(in) :: BLaurent
      real(kind(1d0)) , dimension(:,:,-4:) , intent(in) :: BijLaurent
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
    end subroutine CalcI1nloddim
  end interface
!
  dsigmaB  = 0
  dsigmaV  = 0
  dsigmaVV = 0
!
! if the weight is zero or both of the contributions are
! not needed return back:
  if ((weightPS.eq.0).or. &
      ((.not.flg_LO).and. &
       (.not.flg_NLO_V).and. &
       (.not.flg_NNLO_VV))) then
    return
  end if
!
! Calculation of the flux factor:
  fluxfact = calcFluxFact(p)
! Dealing only with tree-level amplitudes hence no scale
! dependence occurs if couplings are factored out.
! If only the Born is needed the Born SME is calculated only:
  if (flg_LO.and..not.(flg_NLO_V.or.flg_NNLO_VV)) then
    call calcB(p,smeB)
! If the virtual is needed the corresponding I operator needs
! the color-correlated Born SMEs, hence the Born SME is obtained from
! the color-correlated one:
  elseif (flg_NLO_V.and..not.flg_NNLO_VV) then
!    call calcBij(p,Bij)
! Uncomment the following four lines if you want to use the d-dimensional
! Born:
    call CalcBijddim(p,Bij,BijLaurent)
    do i=-4,2
      call CastSMEijToSME(p,BijLaurent(:,:,i),BornLaurent(i))
    end do
    call CastSMEijToSME(p,Bij,smeB)
  elseif (flg_NNLO_VV) then
!    call CalcBijkl(p,Bijkl)
    call CalcBijklddim(p,Bijkl,BijklLaurent)
    do i=-4,2
      call CastSMEijklToSMEij(.false.,p, &
                              BijklLaurent(:,:,:,:,i),BijLaurent(:,:,i))
      call CastSMEijToSME(p,BijLaurent(:,:,i),BornLaurent(i))
    enddo
    smeB = BornLaurent(0)
    Bij  = Bijlaurent(:,:,0)
  end if
! Storing original xi values:
!  print *,"xir,xif: ",xir,xif
  call StoreXis_scales
! Run through all the scales:
  do iscale=1,nscales
! Changing the scales and recalculating the couplings accordingly:
    call ChangeXis_scales(iscale)
! Setting up the new muR and muF:
    call calcmyscales(p)
! Recalculating the couplings:
    call calc_couplings
! Calling the virtual routines (if needed):
    if (flg_NLO_V.and..not.flg_NNLO_VV) then
!      call CalcV(p,smeV)
!      call CalcI1nlo(p,smeB,Bij,I_NLO)
! Virtual present up to O(ep^2):
      call CalcVddim(p,smeV)
      call CalcI1nloddim(p,BornLaurent,BijLaurent,I_NLO)
    elseif (flg_NNLO_VV) then
!      call CalcVij(p,Vij,VijLaurent)
      call CalcVijddim_new(p,Vij,VijLaurent)
      call CastSMEijToSME(p,Vij,smeV)
!      call CalcI1nlo(p,smeB,Bij,I_NLO,Laurent1)
      call CalcInnlo(p,                  &
                     smeB,BornLaurent,   &
                     Bij,BijLaurent,     &
                     Bijkl,BijklLaurent, &
                     smeV,VijLaurent,    &
                     I_NLO,I_NNLO,I_NLOLaurent)
      call CalcVV(p,smeVV)
    end if
! Obtaining the Born contribution:
    if (flg_LO) then
      dsigmaB(iscale) = smeB &
                      * weightPS &
                      * fluxfact &
                      * alphas**border_as &
                      * alphaEM**border_aEM
    end if
! Obtaining the Virtual + I1:
    if (flg_NLO_V) then
      dsigmaV(iscale) = (smeV + I_NLO) &
                      * weightPS &
                      * fluxfact &
                      * alphas**(border_as+1) &
                      * alphaEM**border_aEM
    end if
! Obtaining the contribution from I operator(s) at NNLO:
    if (flg_NNLO_VV) then
      dsigmaVV(iscale) = (smeVV + I_NNLO) &
                       * weightPS &
                       * fluxfact &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
    end if
  end do
! Restoring the old xi values:
  call RestoreXis_scales
!
!
end subroutine sigma_B_V_Is
