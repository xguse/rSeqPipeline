from rSeq.utils.sitRep import start_sitrep
start_sitrep('/home/augustine/tmp/')
print "this should be stdOut!\nBut also this!\n"
raise Exception('This should be stdErr!')