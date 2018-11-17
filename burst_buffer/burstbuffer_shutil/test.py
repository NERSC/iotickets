import shutil
file1=sys.argv[1]
file2=sys.argv[2]
try:
  shutil.copy2(file1,file2)
except Exception as e:
  print ('error in copy2:',e)

try:
  shutil.copy(file1,file2)
except Exception as e:
 print ('error in copy:',e)
