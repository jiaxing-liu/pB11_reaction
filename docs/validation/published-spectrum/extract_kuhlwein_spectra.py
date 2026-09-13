import xml.etree.ElementTree as ET,re,csv,json,hashlib,sys
from pathlib import Path
root=ET.parse(sys.argv[1]).getroot()
curves={}
for n in root.iter():
 if n.get('stroke-width')!='1.5' or not n.get('transform','').endswith('316.262348, 666.539298)'):continue
 nums=[float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?',n.get('d',''))]
 if len(nums)<4:continue
 pts=list(zip(nums[::2],nums[1::2]));my=sum(y for x,y in pts)/len(pts)
 panel=0 if my>490 else 1 if my>260 else 2
 color='data' if n.get('stroke')=='rgb(0%, 0%, 0%)' else 'model'
 curves.setdefault((panel,color),[]).extend(pts)
with open(sys.argv[2],'w') as f:
 w=csv.writer(f,lineterminator="\n");w.writerow(['panel','curve','energy_keV','relative_height'])
 for (panel,color),pts in sorted(curves.items()):
  y0,y1=[(493.602868,696.765556),(266.670318,469.827588),(34.828835,241.719279)][panel]
  byx={x:y for x,y in pts}
  for x,y in sorted(byx.items()):w.writerow([panel,color,6000*(x-43.220504)/(273.734516-43.220504),(y-y0)/(y1-y0)])
print({str(k):len(set(x for x,y in v)) for k,v in curves.items()})
