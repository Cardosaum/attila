f = r"cdr3,quantity,length,MW,AV,IP,flex,gravy,SSF_Helix,SSF_Turn,SSF_Sheet,nºA,nºC,nºD,nºE,nºF,nºG,nºH,nºI,nºK,nºL,nºM,nºN,nºP,nºQ,nºR,nºS,nºT,nºV,nºW,nºY,%A,%C,%D,%E,%F,%G,%H,%I,%K,%L,%M,%N,%P,%Q,%R,%S,%T,%V,%W,%Y"

for r in range(f.count(',')):
  f = f.replace(',',';')
print(f)
