import dict_utils
import sys
import os

print "Reading dict file"
p = dict_utils.flipperDict()
p.read_from_file(sys.argv[1])

nSims=p['nSims']
for iii in range(nSims):
    os.system('python get_spectra.py %s %03d'%(sys.argv[1],iii))
    os.system('python process_spectra.py %s %03d'%(sys.argv[1],iii))



