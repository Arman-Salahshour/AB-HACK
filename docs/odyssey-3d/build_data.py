import re,html,json
s=open("../final-ranking.html").read()
s=re.sub(r'<style.*?</style>|<script.*?</script>|<svg.*?</svg>','',s,flags=re.S)
s=re.sub(r'data:[^"\')]{200,}','DATA',s)
t=re.sub(r'<tr[^>]*>','\n',s); t=re.sub(r'<[^>]+>',' ',t); t=html.unescape(re.sub(r'[ \t]+',' ',t))
rows=[]
for l in t.split('\n'):
    m=re.match(r'\s*(\d+) (\S+) (?:(twin)|(optimisation|de novo)) (?:(?:optimisation) )?(\S+) (.*)',l)
    if not m or int(m.group(1))>30: continue
    r,name,twin,track,cdr,rest=m.groups()
    track=track or 'optimisation'
    d=dict(rank=int(r),name=name,twin=bool(twin),track=track,cdr=cdr)
    if track=='optimisation':
        n=rest.strip().split()[-6:]
        d.update(prob=float(n[0]),kd=float(n[1]),ipae=float(n[2]),pose=float(n[3]),ph=float(n[4]),gap=float(n[5]))
    else:
        m2=re.search(r'ipTM ([\d.]+), interface PAE (\d+)',rest)
        d.update(prob=0,iptm=float(m2.group(1)),ipae=float(m2.group(2)))
    rows.append(d)
print(len(rows))
json.dump(rows,open("data.json","w"))
