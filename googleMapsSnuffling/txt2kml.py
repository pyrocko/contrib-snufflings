from pyrocko import guts
import re

f = open('Plates.txt', 'r')

sections = []
headers = []
boundary = []
for l in f.readlines():
    if l.startswith('*'):
        sections.append(boundary)
        boundary = []
    
    elif re.match('[A-Z]', l):
        headers.append(l)

    else:
        try:
            boundary.append(map(float, l.split(',') ))
        except ValueError:
            pass
f.close()


datastr = ''
for i, bound in enumerate(sections):
    header = headers[i]
    datastr += '<Placemark>\n<name>%s</name>\n<LineString>\n<coordinates>\n'%\
            " ".join(header.split()[1::])

    for segment in bound:
        print tuple(segment)+tuple([100])
        datastr += "%s, %s, %s\n"%(tuple(segment)+tuple([100.]))
    datastr += '</coordinates>\n</LineString>\n</Placemark>\n'


filestr = '''
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.0"> 
<Document>
%s
</Document>
</kml>
'''%datastr

outfn = 'plates_withnames.kml'
fkml = file(outfn, 'w')
fkml.write(filestr)
fkml.close()
