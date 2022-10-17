
# Makefile to create all PoWR exe-files (.exe and .exe.opt) from the .f files in the same directory
# H. Todt 09.10.2022 htodt@astro.physik.uni-potsdam.de

linkdate=$(shell date)
linkuser=$(shell whoami)
linkhost=$(shell hostname)

FC = ifort -i8 -r8
FFLAGS     = -assume byterecl -save -extend-source  -O1 -fpe0 -traceback -mcmodel medium -g -fpconstant -fp-model strict
FFLAGS_OPT = -assume byterecl -save -extend-source  -O2 -fpe0 -traceback -mcmodel medium -g -axCORE-AVX2 -ipo -unroll -fpconstant -fp-model strict

MKLPATH    = ${MKLROOT}/lib/intel64
MKLINCLUDE = ${MKLROOT}/include

LINKER_OPTIONS = -L${MKLPATH} -I${MKLINCLUDE} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -mcmodel medium 
LINKER_DYNAMIC = -shared-intel 
# LINKER_STATIC  = -L${MKLPATH} -I${MKLINCLUDE} -Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_sequential.a ${MKLPATH}/libmkl_core.a -Wl,--end-group -lpthread 

ADAPTEROBJECTS  =  adapop.o  adapter.o  adatrans.o  addhistentry.o  chrinstr.o  clock.o  closms.o  cmsstore.o  count.o  datom.o  decadp.o  fedat.o  findcharge.o  idx.o  install.o  isrcheq.o  jsymset.o  lengthms.o  lipo.o  openms.o  priadp.o  readms.o  remark.o  rmodadp.o  sargc.o  sargp.o  sargrest.o  sargv.o  second.o  stamp.o  storage.o  writms.o 
COLIOBJECTS     =  addhistentry.o  addopa.o  backjc.o  bnue.o  calc_s.o  check_cont.o  check_lines.o  cldiffus.o  clloade.o  clock.o  clopene.o  closms.o  clsavee.o  clsavewc.o  cmfcoop.o  cmffeop.o  cmsstore.o  coli.o  colihist.o  colimo.o  colimop.o  coli_setzero.o  coliwm.o  coop.o  count.o  datom.o  dbnuedt.o  deccoli.o  difdtdr.o  drlevel.o  fecheck.o  fedat.o  findcharge.o  folr.o  frequint.o  frequnorm.o  gauntff.o  genwp1.o  gethistentry.o  horner.o  idx.o  install.o  inv.o  invtri.o  isamax.o  isrcheq.o  isrchfgt.o  jsymset.o  ksigma.o  liop.o   my_isamax.o  opaross.o  openexms.o  openms.o  owninv.o  photocs.o  photon3.o  plotalpha.o  plotanf.o  plotcon.o  plotcons.o   polyfit.o  popmin_nulling.o  preline.o  prepk.o  prep_ppp.o  pri_epsg.o  readms.o  remark.o  rmodcoli.o  rsort.o  sargc.o  sargp.o  sargv.o  second.o  seqlinecl.o  set_momzero.o  setup_ff.o  sfit.o  sfitplo.o  shift.o  shortray.o  splinpo_fast.o   stamp.o  storage.o  store_ff.o  tradfun.o  vmf.o  wmodcoli.o  writms.o 
COMOOBJECTS     =  addhistentry.o  bfcross.o  bnue.o  clock.o  closms.o  cmsstore.o  como.o  coop.o  count.o  datom.o  dbnuedt.o  decnot.o  decomo.o  difdtdr.o  diffus.o  fedat.o  findcharge.o  folr.o  gauntff.o  gethistentry.o   idx.o  install.o  invtri.o  isrcheq.o  jsymset.o  ksigma.o  momo.o  opaross.o  openms.o  photocs.o  photon3.o  plohtot.o  plotanf.o  plotanfs.o  plotcon.o  plotcons.o  plotopa.o  plotrtau1.o   popmin_nulling.o  prihtot.o  primint.o  priopa.o  readms.o  remark.o  remoco.o  sargc.o  sargp.o  sargv.o  second.o  stamp.o  storage.o  tradfun.o  writms.o 
EXTRAPOBJECTS   =  addhistentry.o  aitken.o  change.o  clock.o  closms.o  cmsstore.o  count.o  datom.o  decnot.o  extrap.o  fedat.o  findcharge.o  gauntff.o   idx.o  inhibit.o  isrcheq.o  jsymset.o  ng3.o  ng4.o  openms.o  pricorr.o  priex.o  readms.o  remark.o  sargc.o  sargp.o  sargv.o  second.o  stamp.o  storage.o  writms.o 
FORMALOBJECTS   =  backjc.o  bandwidth.o  bnue.o  clock.o  closms.o  cmffeop.o  cmfset_formal.o  cmsstore.o  convolgauss_flex.o  convolopafe.o  coop.o  copy_secondmodel.o  count.o  datom.o  dbnuedt.o  decform.o  diffus.o  drtrans.o  elimin.o  equal.o  fecheck.o  fedat.o  fierfc.o  filterfunctions.o  findcharge.o  findind.o  findldr.o  folr.o  formal.o  formcmf.o  formosa.o  gauntff.o  genw0.o   gmalu.o  horner.o  idx.o  insert_line.o  install.o  inv.o  invtri.o  isamax.o  isrcheq.o  isrchfge.o  isrchfgt.o  isrchfle.o  isrchflt.o  jsymset.o  kholtsmark.o  ksigma.o  limbdark_output.o  limbdark_prep.o  linstark.o  liop.o  lipo.o  ltepop.o  macroclump.o   manipop.o  mdmv.o  mdv.o  merge_rgrid.o  moment0.o  moment1.o  moment2.o  msub.o  multiple.o  multspli.o  mvmd.o  mvv.o  my_isamax.o  newpol2.o  nowind.o  obsfram.o  openms.o  owninv.o  phiholtsmark.o  photocs.o  photon3.o  plotanf.o  plotanfs.o  plotcon.o  plotcons.o  plot_secondmodel_grid.o  plotvdop.o  plot_windrot_grid.o  polyfit.o  popmin_nulling.o  preform.o  prepmacroclump.o  prepray.o  pridwl.o  priopal.o  pri_par.o  pripro.o  quadstark.o  read_h_starkdata.o  read_linecard_parameters.o  readms.o  remark.o  rescale_secmod.o  rotation_prep.o  sargc.o  sargp.o  sargrest.o  sargv.o  sdot.o  second.o  secondmodel_define.o  secondmodel_prep.o  set_pop_zero.o  setup.o  sfit.o  sfitplo.o  shift.o  sofbet.o  splinpo.o  splinpo_fast.o  stamp.o  starkbroad.o  starkdamp_hei.o  starkheii.o  starkheiiprep.o  starkhi.o  stark_hi_lemke.o  starkholtsmark.o  starkprof.o  starkvoigt.o  storage.o  taucmu.o  tradfun.o  tradwl.o  transform_rgrid.o  traplo.o  trbk.o  vadd.o  vcse1f.o  vdop_struct.o  vmalv.o  vmf.o  voigth.o  writms.o  zoneint.o 
MODIFYOBJECTS   =  addhistentry.o  change.o  clock.o  closms.o  cmsstore.o  count.o  datom.o  fedat.o  findcharge.o  gethistentry.o  idx.o  inhibit.o  intepo.o  isrcheq.o  jsymset.o  modify.o  openms.o  primo.o  readms.o  remark.o  sargc.o  sargp.o  sargv.o  second.o  stamp.o  storage.o  writms.o 
MSINFOOBJECTS   =  idx.o  msinfo.o  storage.o
NEWDATOMOBJECTS =  findelement.o  idx.o  jsymset.o  newdatom.o  newdatomion.o  sargc.o  sargp.o  sargv.o 
NEWFORMAL_CARDSOBJECTS =  clock.o  closms.o  cmsstore.o  count.o  datom.o  fedat.o  findcharge.o  idx.o  install.o  isrcheq.o  jsymset.o  newformal_cards.o  openms.o  readms.o  remark.o  sargc.o  sargp.o  sargrest.o  sargv.o  storage.o 
NJNOBJECTS     =  addhistentry.o  closms.o  cmsstore.o  idx.o  jsymset.o  njn.o  openms.o  readms.o  storage.o  writms.o 
STEALOBJECTS   =  addhistentry.o  addopa.o  adjgamma.o  aitken.o  backuppopnum.o  bfcross.o  bnue.o  brmv.o  brnorm2.o  brtpup.o  brvdivs.o  brvm.o  brvvdy.o  calcgammaradmean.o  calcmassfromgeff.o  calcwc.o  cbbfe.o  cbbh.o  cbbhe.o  cbbmore.o  cbbn.o  ccore.o  change.o  clock.o  closms.o  clump_struct.o  cmsstore.o  cofreq.o  colli.o  coma.o  coop.o  coopfrq.o  count.o  datom.o  dbnuedt.o  dcoop.o  decnot.o  decste.o  decvelpar.o  delpla.o  deltagr.o  deltagrthin.o  deriv.o  dliop.o  dmopen.o  drdat.o  ensuretaumax.o  erf.o  expint1exp.o  extrap.o  fecheck.o  fedat.o  feop_steal.o  filterfunctions.o  findcharge.o  flag_zerorates.o  flgrid.o  funsch.o  gauntff.o  geomesh.o  gethistentry.o  gradiff.o  hysthdruku.o  hystruku.o  idx.o  inhibit.o  initfcorr.o  initvel.o  initvelbetapar.o  install.o  interpolatepopnum.o  interpolatetemp.o  inv.o  isamax.o  ismax.o  isrcheq.o  isrchfge.o  isrchfgt.o  isrchfle.o  isrchflt.o  isrchne.o  jlderiv.o  jsymset.o  ksigma.o  lcore.o  lengthms.o  linpop.o  linsol.o  linsol_split.o  liop.o  lipo.o  load_ff.o  loadwc.o  ltepop.o  mgoetz.o  mlanger.o  my_isamax.o  nextjob.o  ng3.o  ng4.o  nltepop.o  opaross.o  openexms.o  openms.o  overlap.o  owninv.o  pgrid.o  photocs.o  photon3.o  plocc.o  plotacc.o  plotalpha.o  plotanf.o  plotanfs.o  plotapp.o  plotcon.o  plotcons.o  plotdep.o  plotflu.o  plotgamma.o  plothsum.o  plotjline.o  plotjnue.o  plotpop.o  plotsigmafe.o  plott.o  plotunlu.o  plotv.o  plotvgrad.o  popzero.o  preline.o  prepkubat.o  preplotapp.o  pri1rat.o  pricc.o  pricolr.o  pricomp.o  pricorr.o  pridat.o  priex.o  priexpo.o  priflux.o  prigahist.o  prihist.o  prijost.o  prilc.o  primat.o  primod.o  printmodelsummary.o  pripop.o  prirat.o  pritau.o  priunlu.o  radio.o  radnet.o  readms.o  redcor.o  regridoldpop.o  regula.o  remark.o  remarkf.o  remost.o  rgrid.o  rsort.o  sargc.o  sargp.o  sargv.o  sdot.o  second.o  seqlinecl.o  setxjc.o  setxjfine.o  setxjl.o  shift.o  shiftreal.o  shiftstring.o  smach.o  splinpo.o  splinpox.o  stamp.o  steal.o  sthist.o  storage.o  tauscal.o  tcolor.o  tdiffus.o  tempcorr.o  tempcorr_expdamp.o  tempcorr_fluxerr.o  tradfun.o  trbk.o  velobeta.o  velthin.o  vmf.o  vsub.o  writms.o  wrvel.o  xextinc.o  xrudi.o  zanstra.o 
WRCONTOBJECTS  =  addhistentry.o  bnue.o  clock.o  closms.o  cmsstore.o  coop.o  count.o  datom.o  dbnuedt.o  decon.o  delpla.o  difdtdr.o  diffus.o  elimin.o  equal.o  fedat.o  filterfunctions.o  findcharge.o  gauntff.o  gethistentry.o   horner.o  idx.o  install.o  inv.o  isamax.o  isrcheq.o  jsymset.o  ksigma.o  lipo.o  mdmv.o  mdv.o  moment0.o  moment1.o  moment2.o  msub.o  mvmd.o  mvv.o  my_isamax.o  opaross.o  openms.o  openmsr.o  owninv.o  photocs.o  photon3.o  plotanf.o  plotcon.o  plotcons.o  polyfit.o  popmin_nulling.o  pricolr.o  priflux.o  priint.o  prijost.o  priopa.o  readms.o  regula.o  remark.o  rmodcon.o  sargc.o  sargp.o  sargv.o  second.o  setup.o  sfit.o  sfitplo.o  splinpo.o  stamp.o  storage.o  tcolor.o  tradfun.o  trbk.o  tremain.o  vadd.o  vmf.o  wrcont.o  writms.o  xextinc.o  xrudi.o  zanstra.o 
WRSTARTOBJECTS =  addhistentry.o  bfcross.o  bnue.o  clock.o  closms.o  clump_struct.o  cmsstore.o  coop.o  coopfrq.o  count.o  datom.o  dbnuedt.o  decfreq.o  decstar.o  decvelpar.o  deltagr.o  deltagrthin.o  dtdr.o  fedat.o  fgrid.o  findcharge.o  gauntff.o  geomesh.o  gradiff.o  grey.o  hysthdruku.o  hystruku.o  idx.o  initvel.o  initvelbetapar.o  install.o  isrcheq.o  isrchfgt.o  isrchne.o  jstart.o  jsymset.o  ksigma.o  lipo.o  lowercase.o  ltepop.o  mgoetz.o  mlanger.o  my_clock.o  my_date.o  opagrey.o  opaross.o  openms.o  pgrid.o  photocs.o  photon3.o  plotanf.o  plotanfs.o  plotcon.o  plotcons.o  plotv.o  plotvgrad.o  prep_gammarad.o  pricomp.o  primod.o  priparam.o  prixdat.o  readms.o  readoldt.o  regula.o  remark.o  rgrid.o  ruku.o  sargc.o  sargp.o  sargv.o  second.o  sequin.o  shiftreal.o  shiftstring.o  splinpo.o  splinpox.o  stamp.o  storage.o  tabread.o  tauscal.o  trbk.o  velobeta.o  velthin.o  writms.o  wrstart.o  wrvel.o 

.PHONY : all all_exe all_opt clean

all : all_exe all_opt
all_exe : adapter.exe coli.exe como.exe extrap.exe formal.exe modify.exe msinfo.exe newdatom.exe newformal_cards.exe njn.exe steal.exe wrcont.exe wrstart.exe
all_opt : adapter.exe.opt coli.exe.opt como.exe.opt extrap.exe.opt formal.exe.opt modify.exe.opt msinfo.exe.opt newdatom.exe.opt newformal_cards.exe.opt njn.exe.opt steal.exe.opt wrcont.exe.opt wrstart.exe.opt

adapter.exe : $(ADAPTEROBJECTS) mainadapter.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

adapter.exe.opt : mainadapter.for $(ADAPTEROBJECTS) 
	$(shell cat mainadapter.for $(patsubst %.o, %.f,$(ADAPTEROBJECTS)) >> all.for )
	$(shell mv all.for mainadapter.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainadapter.for


coli.exe : $(COLIOBJECTS) maincoli.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

coli.exe.opt : maincoli.for $(COLIOBJECTS) 
	$(shell cat maincoli.for $(patsubst %.o, %.f,$(COLIOBJECTS)) >> all.for )
	$(shell mv all.for maincoli.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) maincoli.for


como.exe : $(COMOOBJECTS) maincomo.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

como.exe.opt : maincomo.for $(COMOOBJECTS) 
	$(shell cat maincomo.for $(patsubst %.o, %.f,$(COMOOBJECTS)) >> all.for )
	$(shell mv all.for maincomo.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) maincomo.for


extrap.exe : $(EXTRAPOBJECTS) mainextrap.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

extrap.exe.opt : mainextrap.for $(EXTRAPOBJECTS) 
	$(shell cat mainextrap.for $(patsubst %.o, %.f,$(EXTRAPOBJECTS)) >> all.for )
	$(shell mv all.for mainextrap.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainextrap.for


formal.exe : $(FORMALOBJECTS) mainformal.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

formal.exe.opt : mainformal.for $(FORMALOBJECTS) 
	$(shell cat mainformal.for $(patsubst %.o, %.f,$(FORMALOBJECTS)) >> all.for )
	$(shell mv all.for mainformal.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainformal.for


modify.exe : $(MODIFYOBJECTS) mainmodify.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

modify.exe.opt : mainmodify.for $(MODIFYOBJECTS) 
	$(shell cat mainmodify.for $(patsubst %.o, %.f,$(MODIFYOBJECTS)) >> all.for )
	$(shell mv all.for mainmodify.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainmodify.for


msinfo.exe : $(MSINFOOBJECTS) mainmsinfo.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

msinfo.exe.opt : mainmsinfo.for $(MSINFOOBJECTS) 
	$(shell cat mainmsinfo.for $(patsubst %.o, %.f,$(MSINFOOBJECTS)) >> all.for )
	$(shell mv all.for mainmsinfo.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainmsinfo.for


newdatom.exe : $(NEWDATOMOBJECTS) mainnewdatom.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

newdatom.exe.opt : mainnewdatom.for $(NEWDATOMOBJECTS) 
	$(shell cat mainnewdatom.for $(patsubst %.o, %.f,$(NEWDATOMOBJECTS)) >> all.for )
	$(shell mv all.for mainnewdatom.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainnewdatom.for


newformal_cards.exe : $(NEWFORMAL_CARDSOBJECTS) mainnewformal_cards.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

newformal_cards.exe.opt : mainnewformal_cards.for $(NEWFORMAL_CARDSOBJECTS) 
	$(shell cat mainnewformal_cards.for $(patsubst %.o, %.f,$(NEWFORMAL_CARDSOBJECTS)) >> all.for )
	$(shell mv all.for mainnewformal_cards.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainnewformal_cards.for


njn.exe : $(NJNOBJECTS) mainnjn.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

njn.exe.opt : mainnjn.for $(NJNOBJECTS) 
	$(shell cat mainnjn.for $(patsubst %.o, %.f,$(NJNOBJECTS)) >> all.for )
	$(shell mv all.for mainnjn.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainnjn.for

steal.exe : $(STEALOBJECTS) mainsteal.for
	$(FC) -o $(@) $(FFLAGS)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)


steal.exe.opt : mainsteal.for $(STEALOBJECTS) 
	$(shell cat mainsteal.for $(patsubst %.o, %.f,$(STEALOBJECTS)) >> all.for )
	$(shell mv all.for mainsteal.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainsteal.for


wrcont.exe : $(WRCONTOBJECTS) mainwrcont.for
	$(FC) -o $(@) $(FFLAGS)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

wrcont.exe.opt : mainwrcont.for $(WRCONTOBJECTS) 
	$(shell cat mainwrcont.for $(patsubst %.o, %.f,$(WRCONTOBJECTS)) >> all.for )
	$(shell mv all.for mainwrcont.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainwrcont.for


wrstart.exe : $(WRSTARTOBJECTS) mainwrstart.for
	$(FC) -o $(@) $(FFLAGS)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

wrstart.exe.opt : mainwrstart.for $(WRSTARTOBJECTS) 
	$(shell cat mainwrstart.for $(patsubst %.o, %.f,$(WRSTARTOBJECTS)) >> all.for )
	$(shell mv all.for mainwrstart.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) mainwrstart.for

# for each program a file main${progam}.for is required in which ${progam}.f is called
%.for :
	@printf "      PROGRAM MAIN$(patsubst main%.for,%,$(@)) \n"                > main$(patsubst main%.for,%,$(@)).for
	@printf "C***  Provide Link data for possible use in the programm\n"      >> main$(patsubst main%.for,%,$(@)).for
	@printf "      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60\n"      >> main$(patsubst main%.for,%,$(@)).for
	@printf "      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST\n" >> main$(patsubst main%.for,%,$(@)).for
	@printf "      LINK_DATE = '$(linkdate)'\n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "      LINK_USER = '$(linkuser)'\n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "      LINK_HOST = '$(linkhost)'\n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "                               \n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "      CALL $(patsubst main%.for,%,$(@)) \n"                      >> main$(patsubst main%.for,%,$(@)).for
	@printf "      END\n"                                                     >> main$(patsubst main%.for,%,$(@)).for


.f.o :
	$(FC) $(FFLAGS) -c $?

clean :
	rm -f *.o
	rm -f *.for
	rm -f *.exe
	rm -f *.exe.opt
