import sys,os
import commands

def removeNonAscii(s):
    return "".join(
        i for i in s if (
            (ord(i)>=40 and ord(i)<=176) or
            ord(i)==11 or
            ord(i)==12 or
            ord(i)==15 or
            ord(i)==40 or
            ord(i)==39 or
            ord(i)==34
            )
        )

SCRIPTDIR=os.path.dirname(sys.argv[0])

if len(sys.argv)<4:
    print "Usage: %s collectionName host path/to/darshan.log.gz"
    sys.exit(1)

col=sys.argv[1]
host=sys.argv[2]
logfile=sys.argv[3]

cmdstr="%s/darshan-parser-mongo %s 2>/dev/null"%(SCRIPTDIR,logfile)

#ankitb_s3d.x.dns_id2876821_3-7-5469-7598519572177302792_26.darshan.gz

#|tr \"\n\" \" \"|tr -cd '\11\12\15\40-\176'

(status, output)=commands.getstatusoutput(cmdstr)
#print cmdstr,status,output

#print status,output
#if output.find("Error: log file contains out of order rank data.")>=0:
if status!=0:
    content="bad:true"
else:
    content=removeNonAscii(output).replace("\n"," ")

print "db.%s.insert({host:'%s',_id:'%s',%s});"%(col,host,os.path.basename(logfile),content)
