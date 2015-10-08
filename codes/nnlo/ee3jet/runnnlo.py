#!/usr/bin/python
# This script is written to automatize the execution of the NNLO code
# and to process all related tasks:
import os
import sys
import multiprocessing
import time
import math
import random

def GenSeeds():
  global proc
  seedfile = "%s-seeds.dat" %(proc)
  if os.path.isfile(seedfile):
    inp = eval('raw_input("%s is already present, would you like to keep it? (y/n) " %(seedfile))')
    if inp is 'y':
      return
  nseeds = eval(raw_input("How many seeds are required: "))
  digits = eval(raw_input("How many digits are needed:  "))
  # We use the UNIX timestamp for the seed:
  seed = int(time.time())
  random.seed(seed)
  outfile = open(seedfile,'w')
  for iseed in range(nseeds):
    outfile.write("%d\n" %(random.randint(10**(digits-1),10**digits)))
  outfile.close()
  return

# This routine launches the code nnlo  
# and write the standard output to out.iproc:
def worker(iproc):
  global proc
  cmd = "echo '%s\n%d' | %s > out.%d" %(proc,iproc,'./nnlo',iproc)
  print "Process no. %d is started" %(iproc)
#  print cmd
  os.system(cmd)
#  time.sleep(10)
  print "Process no. %d is finished" %(iproc)
  return

# This routine manages the parallel runs:
def runjobs(it,num_proc,njobs):
  jobs = []
  for iproc in range(1,num_proc+1):
# Since we can have multiple iterations we have to assign a
# job id for each of the jobs:
    jobid = num_proc*it + iproc
#    print "it,iproc,jobid: ",it,iproc,jobid
# Since we can have multiple iterations and in the last
# we could have more cores than the number of remaining
# jobs we have to stop when we reach the last job:
    if jobid > njobs:
      break
    p = multiprocessing.Process(target=worker, args=(jobid,))
    jobs.append(p)
    p.start()
# When we started the jobs in iteration it we have to wait
# until of them are finished:
  for p in jobs:
    p.join()
  return

def RunCode():
  global proc
  # We determine the number of cores present: 
  num_cores = multiprocessing.cpu_count()
  print "Number of cores available: ", num_cores
  inp = eval('raw_input("Use this as the number of processes? (Y/N) ")')
  if ('y' in inp)  or ('Y' in inp):
    num_proc = num_cores
  else:
    num_proc = eval(raw_input("The number of processes: "))
  print "\n","The number of processes is: ",num_proc
  njobs = eval(raw_input("The total number of jobs to be run: "))
  # We can have more than one iterations if the number of jobs
  # is greater than the number of cores:
  if njobs <= num_proc:
    itmx = 1
  else:
    itmx = njobs / num_proc + 1
  for it in range(itmx):
    #print "iteration: %d" %(it)
    runjobs(it,num_proc,njobs)
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
#      hista[ibin][3] += histb[ibin][3]**2
      hista[ibin][3] += histb[ibin][2]**2
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

# This routine checks whether we have multiple weights or not, if 
# files are present for multiple weights it returns the number of
# weights otherwise only zero:
def countweights(proc,fname0):
  nwgt = 0
# Construct possible filename tags:
  fname = '%s%s-%04d-W1.dat' %(proc,fname0,1)
  if not os.path.isfile(fname):
    return nwgt
# if we are here we found one of the above file, hence we count down
# how many of them we possess:
  while True:
    fname = '%s%s-%04d-W%d.dat' %(proc,fname0,1,nwgt + 1)
# if any of the above is found quit the loop:
    if not os.path.isfile(fname):
      break
# Otherwise increase nwgt:
    nwgt += 1
# end of while loop...
  print "We found %d weights, multiweight mode is activated..." %(nwgt)
# return the number of weights:
  return nwgt

# This routine reads in the datafiles and tries to collect statistics:
def CollectData(proc):
  proc = "%s-" %(proc)
  fname0 = 'output'
# Before we count down the files we check whether we have multiple 
# weights or no:
  nwgt = countweights(proc,fname0)
  cmulti = ''
# We setup an additional variable for the multiweight analysis:
  if nwgt <> 0:
    iwgt = 1
    cmulti = '-W%d' %(iwgt)
# We try to count the datafiles:
  nfiles = 0
  for ifile in range(1,10000):
    fname = '%s%s-%04d%s.dat' %(proc,fname0,ifile,cmulti)
    if os.path.isfile(fname):
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
      fname = '%s%s-%04d%s.dat' %(proc,fname0,ifile,cmulti)
# We examine all possible files and only consider one if it is 
# present:
      if os.path.isfile(fname):
# We try to open this file:
        if os.path.isfile(fname):
          datafile = open(fname,'r')
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
#        histodata[ihisto][ibin][3] = math.sqrt(histodata[ihisto][ibin][3]/(nfiles-1))
        histodata[ihisto][ibin][3] = math.sqrt(abs(histodata[ihisto][ibin][3]/(nfiles) - histodata[ihisto][ibin][2]**2))
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

# This routine is devoted to collect statistics from different runs:
def CollectStat(proc):
  modswitch = eval('raw_input("Would you like to throw out data? (y/n)")')
  flg_drop = False
  if modswitch == 'y' or modswitch == 'Y':
    flg_drop = True
    threshold = eval(raw_input("Give a threshold factor: "))
  # We loop through all the possible files:
  nfiles = []
  for i in range(1000):
    fname = "%s-stat-%04d.dat" %(proc,i)
    if not os.path.isfile(fname):
      continue
    datafile = open(fname,'r')
    # We have an existing file, if this is the first one we read in the number of
    # contributions calculated:
    if nfiles == []:
      numcont = 0
      conts = []
      while True:
	line = datafile.readline()
	if not line:
	  break
	if '++++' in line or '====' in line or 'xiR' in line:
	  continue
	words = line.split()
	numcont += 1
	cont = [words[0], \
	        float(words[1]), \
		float(words[1])**2, \
		float(words[1])/float(words[2])**2, \
		float(words[1])**2/float(words[2])**2, \
		1/float(words[2])**2]
	conts.append(cont)
	nfiles.append(1)
      # The contributions are already read in we can go to the next file:
      datafile.close()
      continue
    # For all the other files we simply go through all the contributions and add
    # them to the value stored in conts:
    ncont = 0
    while True:
      line = datafile.readline()
      if not line:
	break
      if '++++' in line or '====' in line or 'xiR' in line:
	continue
      words = line.split()
      ncont += 1
      nfiles[ncont-1] += 1
      cont = [words[0], \
	      float(words[1]), \
	      float(words[1])**2, \
	      float(words[1])/float(words[2])**2, \
              float(words[1])**2/float(words[2])**2, \
              1/float(words[2])**2]
      # Check upon correct number of contributions:
      if ncont > numcont: 
	print "Error in %s, more contributions than booked: %d > %d" %(fname,ncont,numcont)
	sys.exit
      # The ordering among contributions should be the same among all the files:
      if cont[0] != conts[ncont-1][0]:
	print "Error in %s, wrong contribution ordering: %s =!= %s" %(fname,cont[0],conts[ncont-1][0])
	sys.exit()
      # If the drop feature is active we calculate the running average from previous runs and
      # compare it with the current file, if the ratio is bigger than the threshold we simply drop 
      # the new file:
      if flg_drop:
	if abs(cont[1]) > threshold*abs(conts[ncont-1][1])/(nfiles[ncont-1] - 1):
	  print "%dth file is dropped for contribution no. %d" %(nfiles[ncont-1],ncont)
          # If current file is dropped we have to decrease the file counter by one unit
	  # for the current contribution:
	  nfiles[ncont-1] -= 1;
	else:
          # Adding contribution of current file:
          conts[ncont-1][1] += cont[1]
          conts[ncont-1][2] += cont[2]
          conts[ncont-1][3] += cont[3]
          conts[ncont-1][4] += cont[4]
          conts[ncont-1][5] += cont[5]
      else:
        # Adding contribution of current file:
        conts[ncont-1][1] += cont[1]
        conts[ncont-1][2] += cont[2]
        conts[ncont-1][3] += cont[3]
        conts[ncont-1][4] += cont[4]
        conts[ncont-1][5] += cont[5]
    datafile.close()
  # Manipulating acquired statistics:
  print "We read in %d files" %(max(nfiles))
  for icont,cont in enumerate(conts):
  # Cross section:
    conts[icont][1] = cont[1] / nfiles[icont]
  # Uncertainty:
    conts[icont][2] = math.sqrt(abs(cont[2] / nfiles[icont] - conts[icont][1]**2) / (nfiles[icont] - 1))
  # chi2:
    conts[icont][3] = (cont[4] - 2*cont[3]*conts[icont][1] + cont[5]*conts[icont][1]**2) / (nfiles[icont] - 1)
  # Print out statistics:
  print "The final statistics: "
  for icont,cont in enumerate(conts):
    print "%s %f +- %f , chi2 = %f" %(cont[0].ljust(6),cont[1],cont[2],cont[3])
  return

if __name__ == '__main__':
  proc = eval('raw_input("The name of the process: ")')
  # List of possible actions:
  modswitch = eval('raw_input("Choose: Full run (F/f), generate random seeds (G/g), collect statistics (S/s), collect histos (H/h) ")')
  if modswitch is 'G' or modswitch is 'g':
    GenSeeds()
  if modswitch is 'F' or modswitch is 'f':
    # there is no input card we quit:
    inpcard = "%s-nnlo.input" %(proc)
    if not os.path.isfile(inpcard):
      print "no input card is found..."
      sys.exit()
    # Check upon random seeds: 
    seedfile = "%s-seeds.dat" %(proc)
    if not os.path.isfile(seedfile):
      print "No file is present with random seeds, generating..."
      GenSeeds()
    # We measure the execution time:
    time_start = time.time()
    # Call the Code runner:
    RunCode()
    time_finish = time.time()
    dt = time_finish - time_start
    print "The elapsed time is: %dh %dmin %ds" %(int(dt)/3600,(int(dt) % 3600) / 60,int(dt) % 60)
  if modswitch is 'H' or modswitch is 'h':
    CollectData(proc)
  if modswitch is 'S' or modswitch is 's':
    CollectStat(proc)
