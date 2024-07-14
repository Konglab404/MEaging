import sys
import gzip

filename = str(sys.argv[1])
print(filename)

f = gzip.open(filename)

f1 = open("1.bis",'w')
f2 = open("2.bis",'w')
f3 = open("3.bis",'w')
f4 = open("4.bis",'w')
f5 = open("5.bis",'w')
f6 = open("6.bis",'w')
f7 = open("7.bis",'w')
f8 = open("8.bis",'w')
f9 = open("9.bis",'w')
f10 = open("10.bis",'w')
f11 = open("11.bis",'w')
f12 = open("12.bis",'w')
f13 = open("13.bis",'w')
f14 = open("14.bis",'w')
f15 = open("15.bis",'w')
f16 = open("16.bis",'w')
f17 = open("17.bis",'w')
f18 = open("18.bis",'w')
f19 = open("19.bis",'w')
f20 = open("20.bis",'w')
f21 = open("21.bis",'w')
f22 = open("22.bis",'w')

line = f.readline()

while line:
  if line.find('chr1\t') != -1:
    f1.writelines(line.replace('chr1\t','1\t'))
  elif line.find('chr2\t') != -1:
    f2.writelines(line.replace('chr2\t','2\t'))
  elif line.find('chr3\t') != -1:
    f3.writelines(line.replace('chr3\t','3\t'))
  elif line.find('chr4\t') != -1:
    f4.writelines(line.replace('chr4\t','4\t'))
  elif line.find('chr5\t') != -1:
    f5.writelines(line.replace('chr5\t','5\t'))
  elif line.find('chr6\t') != -1:
    f6.writelines(line.replace('chr6\t','6\t'))
  elif line.find('chr7\t') != -1:
    f7.writelines(line.replace('chr7\t','7\t'))
  elif line.find('chr8\t') != -1:
    f8.writelines(line.replace('chr8\t','8\t'))
  elif line.find('chr9\t') != -1:
    f9.writelines(line.replace('chr9\t','9\t'))
  elif line.find('chr10\t') != -1:
    f10.writelines(line.replace('chr10\t','10\t'))
  elif line.find('chr11\t') != -1:
    f11.writelines(line.replace('chr11\t','11\t'))
  elif line.find('chr12\t') != -1:
    f12.writelines(line.replace('chr12\t','12\t'))
  elif line.find('chr13\t') != -1:
    f13.writelines(line.replace('chr13\t','13\t'))
  elif line.find('chr14\t') != -1:
    f14.writelines(line.replace('chr14\t','14\t'))
  elif line.find('chr15\t') != -1:
    f15.writelines(line.replace('chr15\t','15\t'))
  elif line.find('chr16\t') != -1:
    f16.writelines(line.replace('chr16\t','16\t'))
  elif line.find('chr17\t') != -1:
    f17.writelines(line.replace('chr17\t','17\t'))
  elif line.find('chr18\t') != -1:
    f18.writelines(line.replace('chr18\t','18\t'))
  elif line.find('chr19\t') != -1:
    f19.writelines(line.replace('chr19\t','19\t'))
  elif line.find('chr20\t') != -1:
    f20.writelines(line.replace('chr20\t','20\t'))
  elif line.find('chr21\t') != -1:
    f21.writelines(line.replace('chr21\t','21\t'))
  elif line.find('chr22\t') != -1:
    f22.writelines(line.replace('chr22\t','22\t'))
  line=f.readline()


f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()
f8.close()
f9.close()
f10.close()
f11.close()
f12.close()
f13.close()
f14.close()
f15.close()
f16.close()
f17.close()
f18.close()
f19.close()
f20.close()
f21.close()
f22.close()
