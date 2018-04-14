from pymprog import *


#fb=open('./data/u2p_score_0319.txt','r')
#nusr=1000
#nmsg=17
fb=open('./data/u2p_score_0412.txt','r')
nusr=10000
nmsg=17

lines=fb.readlines()
Smin=[500 for i in range(nmsg)]
for i in range(10,13):
    Smin[i] = 800
for i in range(13,nmsg):
    Smin[i] = 300
#Smax=

users={}
usersrank={}
items={}
userid={}
usercount=0
itemid={}
itemcount=0
for lin in lines:
    ui,ii,scorei,notknown = lin.split(',')
    if ui not in users:
        users[ui]=set()
        usersrank[ui]=set()
    users[ui].add((ii,float(scorei)))
    usersrank[ui].add((ii,float(notknown)))
    if ii not in items:
        items[ii]=set()
    items[ii].add((ui,float(scorei)))
    if ui not in userid:
        userid[ui]=usercount
        usercount+=1
    if ii not in itemid:
        itemid[ii]=itemcount
        itemcount+=1
scoreui=[[0.0]*nmsg for i in range(nusr)]
#ranki=[[0.0]*17 for i in range(1000) ]
for u,uitems in users.items():
    for i,iscore in uitems:
        scoreui[userid[u]][itemid[i]]=iscore
        #ranki[userid[u]][itemid[i]]=

smax = 0.0
for i in range(nusr):
    for j in range(nmsg):
        if smax < scoreui[i][j]:
            smax = scoreui[i][j]
for i in range(nusr):
    for j in range(nmsg):
        scoreui[i][j] = scoreui[i][j] / smax


begin('wuyan')
irange=iprod(range(len(items)))
urange=iprod(range(len(users)))
x=var('x',iprod(range(len(users)),range(len(items))),kind=bool)

for i in range(len(items)):
    sum(x[j,i] for j in range(len(users)))>=Smin[i]
    #sum(x[j,i] for j in range(len(users)))>=Smin, sum(x[j,i] for j in range(len(users)))<=Smax
for i in range(len(users)):
    sum(x[i,j] for j in range(len(items)))==1

maximize(sum(scoreui[i][j]*x[i,j] for i in range(len(users)) for j in range(len(items))))
solve()


'''
for i in range(len(items)):
    for j in range(len(users)):
        if x[j,i].primal>0:
            print(j,i,x[j,i].primal)
'''

#sensitivity()

asum = [0 for i in range(nmsg)]

fb=open('./output/result1','w')
fb.write('[')
for i in range(len(users)):
    fb.write('[')
    for j in range(len(items)):
        fb.write(str(x[i,j].primal)+" ")
        asum[j] = asum[j] + x[i,j].primal
    fb.write(']\n')
fb.write(']')
fb.close()

for i in range(nmsg):
    print(i, asum[i])

print("total score: "+str(sum(scoreui[i][j]*x[i,j].primal for i in range(len(users)) for j in range(len(items)))/(nusr * nmsg)))

