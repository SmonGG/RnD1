def arraysum(a1, a2):
	sum=0
	i=0
	while i < len(a1) or i < len(a2):
		sum = sum + a1[i]+a2[i]
		i = i + 1
	print sum

def arraytimes(a1, fac):
	sum=0
	i=0
	while i < len(a1):
		a1[i] = a1[i] * fac
		i = i + 1
	print a1

def correction(theo, dist):
	dif = arraysum(theo, dist)
	print dif
	soln = False
	factor = 10
	dist2 = arraytimes(dist, 10)
	corr = 0
	while soln == False or corr < 1000:
		dif2 = arraysum(theo, dist2)
		print dif2
		if dif > dif2:
			dif = dif2
			corr = corr + factor
		else:
			soln = True
		soln = True
	return corr

a1=[100,200,400,800,400,200,100]
a2=[10,20,40,80,40,20,10]
print correction(a1, a2)