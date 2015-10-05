#!/usr/bin/python
# This script is written to help running the analysis (PYTHIA,HERWIG,lhef_analysis)
# in a parallel fashion than combine the various histo files into one:
import os
import sys
import multiprocessing
import time
import math

# This routine launches the code analprog to run on eventfile 
# and write the standard output to out.iproc:
def worker(iproc,nicetag,analprog,eventfile):
  global flg_powhel
  global modeswitch
  global cweight
  global multi_low,multihigh
  if not flg_powhel:
    cmd = "echo '%s' | %s %s > out.%d" %(eventfile,nicetag,analprog,iproc)
  if flg_powhel and modeswitch == 0:
    cmd = "echo '%d\n%s' | %s %s > out.%d" %(modeswitch,eventfile,nicetag,analprog,iproc)
  if flg_powhel and modeswitch > 0 and modeswitch != 3:
    cmd = "echo '%d\n%d\n%s' | %s %s > out.%d" %(modeswitch,cweight,eventfile,nicetag,analprog,iproc)
  if flg_powhel and modeswitch == 3:
    cmd = "echo '%d\n%d\n%d\n%s' | %s %s > out.%d" %(modeswitch,multi_low,multi_high,eventfile,nicetag,analprog,iproc)
  print "Process no. %d is started" %(iproc)
#  print cmd
  os.system(cmd)
#  time.sleep(10)
  print "Process no. %d is finished" %(iproc)
  return

# This routine manages to run the jobs parallel, to work
# we have to give it the code to call (analprog), the event
# file template (evttemplate) and the iteration (it):
def runjobs(it,num_proc,nfiles,nicetag,analprog,evttemplate):
  jobs = []
  for iproc in range(1,num_proc+1):
# Since we can have multiple iterations we have to assign a
# job id for each of the jobs:
    jobid = num_proc*it + iproc
#    print "it,iproc,jobid: ",it,iproc,jobid
# Since we can have multiple iterations and in the last
# we could have more cores than the number of remaining
# jobs we have to stop when we reach the last job:
    if jobid > nfiles:
      break
# We have to set up the event file from the template:
    eventfile = evttemplate + '%04d.lhe' %(jobid)
    p = multiprocessing.Process(target=worker, args=(jobid,nicetag,analprog,eventfile,))
    jobs.append(p)
    p.start()
# When we started the jobs in iteration it we have to wait
# until of them are finished:
  for p in jobs:
    p.join()
  return

def runanalysis(proc,analprog):
# It sometimes happens that the code is not copied to the
# actual folder pt can be found one folder up, this can be
# easily checked:
  if not os.path.isfile(analprog):
    print "%s cannot be found in the folder stepping up..." %(analprog)
    analprog = "../" + analprog
# Check again...
    if not os.path.isfile(analprog):
      print "%s cannot be found..."
      print "Something nasty happens here..."
      sys.exit()
  
# We determine the number of cores in the machine:
  num_cores = multiprocessing.cpu_count()
  print "Number of cores available: ", num_cores

  inp = eval('raw_input("Use this as the number of processes? (Y/N) ")')
  if ('y' in inp)  or ('Y' in inp):
    num_proc = num_cores
  else:
    num_proc = eval(raw_input("The number of processes: "))

  print "\n","The number of processes is: ",num_proc

# We can run the jobs 'nice'-ed:
  nice = eval('raw_input("Would you like the runs as niced (y/n): ")')
  if nice == 'y' or nice == 'Y':
    nicetag = 'nice'
  else:
    nicetag = ''
# We create the template for the event files:
  if proc == '':
# We save the general procid into the proc variable since it is used
# during the analysis routine:
    proc = 'pwg'
    fname0 = 'pwgevents-'
  else:
    fname0 = proc + '-events-'

# We try to count the event files in the current folder:
  print "We count the event files, this can take a while..."
  nfiles = 0
  for ifile in range(1,10000):
    fname = fname0 + '%04d.lhe' %(ifile)
    if os.path.isfile(fname):
      nfiles += 1

# If we found zero event files something definitely went wrong:
  if nfiles == 0:
    print "We failed to find any event file with the %s template" %(fname0)
    sys.exit()

  print "We were able to find %d event files " %(nfiles)

# Now we run the jobs paying attention to the possible 
# iterations arising from the fact that we could have more 
# event files than cores:
  for it in range(nfiles/num_proc + 1):
    runjobs(it,num_proc,nfiles,nicetag,analprog,fname0)
  return

# This routine loads the histogram:
def LoadHisto(datafile,histo):
  while True:
    line = datafile.readline()
# If we hit the end of the file just quit:
    if not line:
      return ''
# We reached a new histogram:
    if '#' in line:
# We store the title of it contained by the line starting with #
      title = line
      break
# We read all the lines until we find an empty one:
  while True:
    line = datafile.readline()
    if not line.strip():  
      break
# It contained real data:
# We have to split the line into segments:
    words = line.split()
# We create a temp list to store the float data:
    templine = []
# We store all the data presented in the actual line in templine
# in float:
    for iword in words:
      templine.append(float(iword))
# The last item should be squared:
    templine[3] = templine[3]**2
# we append the line to the histo:
    histo.append(templine)
  return title

# This routine initializes the histograms:
def InitHistos(datafile,titles,temphistos):
  nhisto = 0
  while True:
# We have to create a temp histogram to hold the actual one:
    temphisto = []
    title = LoadHisto(datafile,temphisto)
    if title == '':
      break
# We append this to temphistos:
    temphistos.append(temphisto)
# We also store the name of the histo:
    titles.append(title)
    nhisto += 1
  print "We loaded %d histograms from the first datafile" %(nhisto)
  return nhisto

# This routine combines histb with hista, the cumulative data is 
# contained by hista:
def CombineHistos(hista,histb):
# We go through all the bins extracting not only the bin content, but
# the ordinal number of the bin as well:
  for ibin, bincontent in enumerate(hista):
# We can combine the bin contents provided by the fact that they
# correspond to the same bin:
    if (bincontent[0] == histb[ibin][0]) and (bincontent[1] == histb[ibin][1]):
      hista[ibin][2] += histb[ibin][2]
      hista[ibin][3] += histb[ibin][3]
    else:
      print "We tried to combine bins, but they mismatch..."
      sys.exit()
  return

# This routine saves the histos with appropriate names to a file:
def SaveHistos(proc,fname0,histos,titles):
# We create the file name:
  fname = "%s-%s.dat" %(proc,fname0)
# if it is present better ask permission:
  if os.path.isfile(fname):
    input = eval("raw_input('%s exists, would you like to overwrite it? (y/n): ' %(fname))")
    if input is 'n' or input is 'N':
      print "bye!"
      sys.exit()
  outfile = open(fname,'w')
# We dump all the histos into the file:
  for ihisto,title in enumerate(titles):
# Title first...
    outfile.write("%s" %(title))
# Then the real data...
    for ibin,bin in enumerate(histos[ihisto]):
      outfile.write("%e   %e   %e   %e\n" %(bin[0],bin[1],bin[2],bin[3]))
# When finished with the histo leave two empty lines:
    outfile.write("\n\n")
  return

# This routine checks weither we have multiple weights or not, if 
# files are present for multiple weights it returns the number of
# weights otherwise only zero:
def countweights(proc,fname0):
  nwgt = 0
# Construct possible filename tags:
  fname1 = '%s-%s-%04d-W1.top' %(proc,fname0,1)
  fname2 = '%s-%04d-%s-W1.top' %(proc,1,fname0)
  fname3 = '%s%s-%04d-W1.top' %(proc,fname0,1)
  fname1 = '%s-%s-%04d-W1.dat' %(proc,fname0,1)
  fname2 = '%s-%04d-%s-W1.dat' %(proc,1,fname0)
  fname3 = '%s%s-%04d-W1.dat' %(proc,fname0,1)
  if not os.path.isfile(fname1) and \
     not os.path.isfile(fname2) and \
     not os.path.isfile(fname3) and \
     not os.path.isfile(fname4) and \
     not os.path.isfile(fname5) and \
     not os.path.isfile(fname6):
    return nwgt
# if we are here we found one of the above file, hence we count down
# how many of them we possess:
  while True:
    fname1 = '%s-%s-%04d-W%d.top' %(proc,fname0,1,nwgt + 1)
    fname2 = '%s-%04d-%s-W%d.top' %(proc,1,fname0,nwgt + 1)
    fname3 = '%s%s-%04d-W%d.top' %(proc,fname0,1,nwgt + 1)
    fname4 = '%s-%s-%04d-W%d.dat' %(proc,fname0,1,nwgt + 1)
    fname5 = '%s-%04d-%s-W%d.dat' %(proc,1,fname0,nwgt + 1)
    fname6 = '%s%s-%04d-W%d.dat' %(proc,fname0,1,nwgt + 1)
# if any of the above is found quit the loop:
    if not os.path.isfile(fname1) and \
       not os.path.isfile(fname2) and \
       not os.path.isfile(fname3) and \
       not os.path.isfile(fname4) and \
       not os.path.isfile(fname5) and \
       not os.path.isfile(fname6):
      break
# Otherwise increase nwgt:
    nwgt += 1
# end of while loop...
  print "We found %d weights, multiweight mode is activated..." %(nwgt)
# return the number of weights:
  return nwgt


# This routine reads in the datafiles and tries to collect statistics:
def collectdata(proc,analtag):
# The datafile format depends on the used analysis:
  if analtag == 'p':
    fname0 = 'POWHEG+PYTHIA-output'
  elif analtag == 'h':
    fname0 = 'POWHEG+HERWIG-output'
  elif analtag == 'l':
    fname0 = 'LHEF_analysis'
  elif analtag == 'n':
    fname0 = 'NLO'
  elif analtag == 'nnlo':
    fname0 = 'output'
# If we did not run the analysis proc can be an empty string:
  if proc == '':
    proc = 'pwg'
# Hack to work even when we really specify the process name:
  if proc is not '':
    proc = "%s-" %(proc)
# Before we count down the files we check weither we have multiple 
# weights or no:
  nwgt = countweights(proc,fname0)
  cmulti = ''
# We setup an additional variable for the multiweight analysis:
  if nwgt <> 0:
    iwgt = 1
    cmulti = '-W%d' %(iwgt)
# We try to count the datafiles:
  nfiles = 0
# We try to count the datafiles two separate ways since the cwhichseed
# tag can be put to two different places depending upon the choice of
# analysis:
  for ifile in range(1,10000):
    fname1 = '%s-%s-%04d%s.top' %(proc,fname0,ifile,cmulti)
    fname2 = '%s-%04d-%s%s.top' %(proc,ifile,fname0,cmulti)
    fname3 = '%s%s-%04d%s.top' %(proc,fname0,ifile,cmulti)
    fname4 = '%s-%s-%04d%s.dat' %(proc,fname0,ifile,cmulti)
    fname5 = '%s-%04d-%s%s.dat' %(proc,ifile,fname0,cmulti)
    fname6 = '%s%s-%04d%s.dat' %(proc,fname0,ifile,cmulti)
    if os.path.isfile(fname1) or \
       os.path.isfile(fname2) or \
       os.path.isfile(fname3) or \
       os.path.isfile(fname4) or \
       os.path.isfile(fname5) or \
       os.path.isfile(fname6):
      nfiles += 1
  print "We found %d datafiles" %(nfiles)
  if nfiles == 0:
    print "We did not find any datafile, cannot proceed..."
    sys.exit()

# loop to treat the case of multiweight case as well:
  while True:
# This list will hold all the histos:
    histodata = []
# This one holds the names for the histos:
    titles = []
# if we have multiple weights we have an extra file tag to apply:
    if nwgt <> 0:
      cmulti = '-W%d' %(iwgt)
# We start to read in the datafiles:
    for ifile in range(1,10000):
      fname1 = '%s-%s-%04d%s.top' %(proc,fname0,ifile,cmulti)
      fname2 = '%s-%04d-%s%s.top' %(proc,ifile,fname0,cmulti)
      fname3 = '%s%s-%04d%s.top' %(proc,fname0,ifile,cmulti)
      fname4 = '%s-%s-%04d%s.dat' %(proc,fname0,ifile,cmulti)
      fname5 = '%s-%04d-%s%s.dat' %(proc,ifile,fname0,cmulti)
      fname6 = '%s%s-%04d%s.dat' %(proc,fname0,ifile,cmulti)
# We examine all possible files and only consider one if it is 
# present:
      if os.path.isfile(fname1) or \
         os.path.isfile(fname2) or \
         os.path.isfile(fname3) or \
         os.path.isfile(fname4) or \
         os.path.isfile(fname5) or \
         os.path.isfile(fname6):
# We try to open one of these files:
        if os.path.isfile(fname1):
          datafile = open(fname1,'r')
        elif os.path.isfile(fname2):
          datafile = open(fname2,'r')
        elif os.path.isfile(fname3):
          datafile = open(fname3,'r')
        elif os.path.isfile(fname4):
          datafile = open(fname4,'r')
        elif os.path.isfile(fname5):
          datafile = open(fname5,'r')
        elif os.path.isfile(fname6):
          datafile = open(fname6,'r')
# When we read in the first file we have to set up the list containing
# all histos:
        if ifile == 1:
          nhisto = InitHistos(datafile,titles,histodata)
# We read in all the histos from the first file, we should add them to
# the histodata array which will hold all the histos, but now it is
# empty:
# If we passed the first file we have to read the histos:
        else:
# we go through all the histos in histodata and combine them with the
# newly read in one:
          for hista in histodata:
    	    histb = []
	    title = LoadHisto(datafile,histb)
	    if title is not '':
	      CombineHistos(hista,histb)
	    else:
	      print "We tried to read in %s or %s, but failed" %(fname1,fname2)
        datafile.close()
# Then we finalize the histograms:
    for ihisto in xrange(len(histodata)):
      for ibin in xrange(len(histodata[ihisto])):
        histodata[ihisto][ibin][2] = histodata[ihisto][ibin][2]/nfiles
        histodata[ihisto][ibin][3] = math.sqrt(histodata[ihisto][ibin][3])/nfiles
# Then, finally, we save the histograms to a file:
    if cmulti == '':
      SaveHistos(proc,fname0,histodata,titles)
    else:
      SaveHistos(proc,fname0 + cmulti,histodata,titles)
# if there is no multiple weight to collect data from or we are at the
# end of the weights simply quit:
    if cmulti == '' or iwgt == nwgt:
      break
# Otherwise go back to the beginning of the loop and increase iwgt by one:
    iwgt += 1
  return

if __name__ == '__main__':
# We ask for the process name or when blank is returned use pwg:
  proc = eval('raw_input("The name of the process: ")')
  if proc == '':
    print "pwg events will be used..."
  anal = eval('raw_input("The type of analysis (Pythia=p,Herwig=h,lhef_analysis=l,NLO=n,NNLO=nnlo): ")')
  if anal == 'p':
    analprog = 'main-PYTHIA-lhef'
  elif anal == 'h':
    analprog = 'main-HERWIG-lhef'
  elif anal == 'l':
    analprog = 'lhef_analysis'
  elif anal == 'n':
    print "This function is only for collecting data..."
  elif anal == 'nnlo':
    analprog = 'nnlo'
  else:
    print "You hit %s, instead of p,h,l or n" %(anal)
    sys.exit()
  flg_anal = eval('raw_input("Would you like to run an analysis? (y/n): ")')
  flg_collect = eval('raw_input("Would you like to collect data? (y/n): ")')
  if flg_anal == 'y' or flg_anal == 'Y':
    powhelmode = eval('raw_input("Would you like to run in PowHel mode (y/n): ")')
    flg_powhel = False
    if powhelmode == 'y' or powhelmode == 'Y':
      flg_powhel = True
    if flg_powhel:
      modeswitch = eval(raw_input("Which mode you want to run the code (0,1,2,3): "))
      if modeswitch > 0 and modeswitch != 3:
	cweight = eval(raw_input("Which weight do you want to use: "))
      if modeswitch == 3:
	multi_low  = eval(raw_input("The first weight which should be considered: "))
	multi_high = eval(raw_input("The last weight which should be considered:  "))
    runanalysis(proc,analprog)
  if flg_collect == 'y' or flg_collect == 'Y':
    collectdata(proc,anal)

