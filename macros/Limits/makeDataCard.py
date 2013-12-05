
import string

shortnames  = ['sig',  'nc', 'qcd', 'gj', 'zj', 'wj', 'wg', 'zg']
theoryunc   = ['1.04', '-',  '1.5', '1.1',   '1.05',  '1.05',  '1.08',   '1.08']
rates       = []

noncolunc  = "1.7"
lumiunc     = "1.026"

systnames   = ['pilup',  'gamma_exp', 'lep_exp', 'jet_exp', 'mc_stat']

###first get the total rate without the signal
print "Now we will compute the expected rate in case of Hypothesis 0 (no signal)"

theRate = 0                                            
with open("systSource.txt") as f:
  content = f.readlines()
  iSample = 0
  for line in content:
    line.strip()
    data = line.split()
    if "ADD" in data[0]: continue
    print data[0],
    print shortnames[iSample]
    theRate = theRate + float(data[1])
    rates.append(data[1])
    iSample =  iSample + 1

print "The expected rate is ",theRate

#now put back the signal rate in the beginning
with open("systSource.txt") as f:
  line = f.readline()
  line.strip()
  data = line.split()
  rates.insert(0, data[1]) 

###prepare the file header                                            
print "imax 1                                                       "
print "jmax *                                                       "
print "kmax *                                                       "
print "------------                                                 "
print "# we have just one channel                                   "
print "bin jet_01                                                   "
print "observation ",theRate,"                                      "
print "                                                             "

###now plug in the rates!
print "bin        ",
for sample in shortnames:
  print "   jet_01",
print ""  

print "process    ",
for name in shortnames:
  print '%9s' % name,
print ""  

iBin = 0
print "process    ",
for name in shortnames:
  print '%9s' % iBin,
  iBin = iBin + 1
print ""  

print "rate       ",
for rate in rates:
  print '%9s' % rate,
print ""  
  
print "------------                                                 "
  
###and now systematics with their nice correlations :) 
print "lumi       ",
print '%7s' % 'lnN',
for sample in shortnames:
  if not "nc" in sample:
    print '%9s' % lumiunc,
  else:
    print '%9s' % "-",
print ""  

print "nc_exp     ",
print '%7s' % 'lnN',
for sample in shortnames:
  if "nc" in sample:
#    print sample,
    print '%9s' % noncolunc,
  else:
    print '%9s' % "-",
print ""  

iSyst = 2
for syst in systnames:
  if not "mc_stat" in syst:
    print '%-11s' % syst,
    print '%7s' % 'lnN',
    
    theRate = 0                                            
    with open("systSource.txt") as f:
      content = f.readlines()
      iSample = 0
      for line in content:
        line.strip()
        data = line.split()
        print '%9s' % data[int(iSyst)],
    print ""
  else:
    with open("systSource.txt") as f:
      content = f.readlines()
      iSample = 0
      for line in content:
        line.strip()
        data = line.split()
        if "nc" in data[0]: 
          iSample = iSample + 1
          continue
        thisSystName = "mc_stat_" + shortnames[iSample] 
        print '%-11s' % thisSystName,
        print '%7s' % 'lnN',
        iName = 0
        for name in shortnames:
#          print name,data[1],rates[iSample]
          if data[1] == rates[iName]:
            print '%9s' % data[int(iSyst)],
          else:
            print '%9s' % "-",
          iName = iName + 1
        print ""
        iSample = iSample + 1
  iSyst = iSyst + 1

iSample = 0
for syst in theoryunc:
  thisSystName = "xsec_" + shortnames[iSample] 
  print '%-11s' % thisSystName,
  print '%7s' % 'lnN',
  iName = 0
  for systval in theoryunc:
    if iName == iSample:
      print '%9s' % systval,
    else:
      print '%9s' % "-",
    iName = iName + 1
  print ""
  iSample = iSample + 1

    
