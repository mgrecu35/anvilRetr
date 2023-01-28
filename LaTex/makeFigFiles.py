# read text from file listOfFiles line by line
# and extract string s from template *{s}*
# and write to file listOfFiles2
# 2019-01-30
f=open("listOfFiles","r").readlines()
#f2=open("listOfFiles2","w")
ifig=1
import os
for line in f:
    line=line.strip()
    # find string s between { and }
    s=line[line.find("{")+1:line.find("}")]
    s2="./Figs/fig%2.2i.png"%ifig
    ifig+=1
    cmd="cp %s %s"%(s,s2)
    print(cmd)
    os.system(cmd)