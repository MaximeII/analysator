import pytools as pt

path_bulk = '/scratch/project_2000203/2D/BFH/reverted_ionosphere_field_boundary/'
path_fig  = '/users/dubartma/analysator/Data_Analysis/test_plot2DFFT/'

pt.plot.plt_2DFFT(filedir=path_bulk,
                  start=200,stop=400,
                  outputdir=path_fig,
                  boxre=[3.0,6.0,15.0,18.0],
                  run='BFH',
                  Bmin=1e-13,Bmax=1e-4,
                  Emin=1e-7,Emax=1e1,
                  nmin=1e3,nmax=1e9,
                  dispersion=['Alfven','FMW','SMW','bulkV'],
                  omegamax=2.8,
                  dr_ref=300000,
                  paper='E')
