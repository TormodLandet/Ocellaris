import os, sys, subprocess

# Run Ocellaris like it was started from the command line
thisdir = os.path.dirname(os.path.abspath(__file__))
inpfile = os.path.join(thisdir, 'lid_driven_cavity_flow.inp')
pyexe = sys.executable
cmd = [pyexe, '-m', 'ocellaris', inpfile]

print 'Running command:\n %s\n' % ' '.join(cmd)
exit(subprocess.call(cmd))
