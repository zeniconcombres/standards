import math
from numpy import log as ln
#####IEC60287 Electric cables – Calculation of the current rating############

####Conductor temperature above ambient
def DT (R,I,Wd,lm1,lm2,n,T1,T2,T3,T4):
	### R : Conductor per unit length AC resistance (Ω/m)
	### I : current (A)
	### Wd : dielectric loss per unit length for the insulation surrounding the conductor (W/m)
	### lm1 : λ1 ratio of losses in the metal sheath to total losses in all conductors in that cable
	### lm2 : λ2 ratio of losses in the armouring to total losses in all conductors in that cable.
	### n : number of load carrying conduttor in the cable
	### T1 : thermal resistance per unit length between one conductor and the sheath (K.m/W)
	### T2 : (If applicable) thermal resistance per unit length of the bedding between sheath and armour(K.m/W)
	### T3 : thermal resistance per unit length of the external serving of the cable (K.m/W)
	### T4 : hermal resistance per unit length between the cable surface and the surrounding medium
	Nj=R*I**2
	P1=(Nj+0.5*Wd)*T1
	P2=(Nj*(1+lm1)+Wd)*n*T2
	P3=(Nj*(1+lm1+lm2)+Wd)*n*(T3+T4)
	DT=P1+P2+P3
	return DT

def CCC(dT, Wd, T1, T2, T3, T4, R_AC, lm1, lm2):
	A=Wd*(0.5*T1+T2+T3+T4)
	B=R_AC*T1
	C=R_AC*(1+lm1)*T2
	D=R_AC*(1+lm1+lm2)*(T3+T4)
	I=((dT-A)/(B+C+D))**0.5
	return I


### Resistance DC - AC
def Rt(R_DC_20,F_a, T):
	Rt=R_DC_20*(1+F_a*(T-20))
	return Rt

def F_skin(Rt, f, ks):
	"""Requires the DC resistance input parameter to be in units Ohm/km
	for this to have the right scale of effect."""
	### F_a : temperature coefficient
	### T : conductor temperature [degree C]
	### ks : Factor in skin effect calculation, check Table 2
	Xs=(((8*math.pi*f)/Rt)*10**(-7)*ks)**0.5
	if Xs <= 2.8:
		F_skin = Xs**4/(192+0.8*(Xs**4))
	elif Xs <= 3.8:
		F_skin=-0.136-0.0177*Xs+0.0563*Xs**2
	else:
		F_skin=0.354*Xs-0.733
	return F_skin

def Proxi_effect(Rt, f, kp,con_OD,s):
	"""Requires the DC resistance input parameter to be in units Ohm/km
	for this to have the right scale of effect."""
	Xp=(((8*math.pi*f)/Rt)*10**(-7)*kp)**0.5 
	# ^this here is square rooted, in paper formula no **0.5
	A=(Xp**4)/(192+0.8*Xp**4)
	B=(con_OD/s)**2
	C=0.312*B
	D=1.18/(((Xp**4)/(192+0.8*Xp**4))+0.27)
	Yp=A*B*(C+D)
	return Yp

def R_AC(Rt,Ys,Yp): ### Conductor AC resistance
	R_AC=Rt*(1+Ys+Yp)
	return R_AC

def R_C(Rt,Ys): ### Conductor AC resistance
	R_C=Rt*(1+Ys)
	return R_C

def t_screen(T,I,R_AC,Wd, T1): ### screen operating temperature
	### T : Conductor maximum operating temperature
	### T1 : thermal resistance between conductor and sheath
	t_s=T-(R_AC*I**2+0.5*Wd)*T1
	return t_s

def Rs(Rs_20,t_screen,F_a):
	Rs=Rs_20*(1+F_a*(t_screen-20))
	return Rs


#### Thermal resistance between conductor and Metallic sheath

def T1 (rt1,t1, con_OD):
	### rt1 : thermal resistivity of insulation
	### t1 : thickness between conductor and sheath
	### con_OD : conductor diameter
	T1=(rt1/(2*math.pi))*ln(1+(2*t1/con_OD))
	return T1

def T2 (rt2,t2, sheath_OD):
	### rt2 : thermal resistivity of beding
	### t2 : insulation thickness
	### sheath_OD : sheath outer diameter
	T2=(rt2/(2*math.pi))*ln(1+(2*t2/sheath_OD))
	return T2

def T3 (rt3,t3, armour_OD):
	### rt3 : thermal resistivity of serving
	### t3 : insulation thickness
	### armour_OD : armour outer diameter
	T2=(rt3/(2*math.pi))*ln(1+(2*t3/armour_OD))
	return T2



def T4_flat_side (rt_s, L, s1, De):
	### rt_s : Soil resistivity
	### L : distance from cable centre to groud surface (mm)
	### s1 : separation of two cable (mm)
	### De : Cable external OD
	u=2*L/De
	A=ln(u+(u**2-1)**0.5)
	B=0.5*ln(1+(2*L/s1)**2)
	T4_side= (rt_s/(2*math.pi))*(A+B)
	return T4_side

def T4_flat_mid(rt_s, L, s1, De):
	u=2*L/De
	A=ln(u+(u**2-1)**0.5)
	B=ln(1+(2*L/s1)**2)
	T4_mid= (rt_s/(2*math.pi))*(A+B)
	return T4_mid

def T4_trefoil(rt_s, L, OD, t3):
	### OD : Cable OD, which different to De
	### Doc : diameter of the imaginary coaxial cylinder which just touchesthe crests of a corrugated sheath
	### t3 : thickness of the serving
	r=OD/2
	Doc = (r/(math.sin(math.pi/6))-r)*2
	De=Doc+2*t3
	u=2*L/De
	T4_tre=(1.5/math.pi)*rt_s*(ln(2*u)-0.63)
	return T4_tre

#### T4 correction

def Gb(Hb,Wb,LG):
	X=min(Hb,Wb)
	Y=max(Hb,Wb)
	A=0.5*(X/Y)*(4/math.pi-X/Y)
	B=ln(1+(Y**2)/(X**2))
	C=ln(X/2)
	rb=math.exp(A*B+C)
	Ub=LG/rb
	Gb= ln(Ub+(Ub**2-1)**0.5)
	return Gb

def T4b(T40,rt_s,rt_b,N,Gb):
	A=T40/rt_s
	B=N/(2*math.pi)*Gb
	T4b=(A-B)*rt_b
	return T4b

def T4s(rt_s,N,Gb):
	T4s=(N/(2*math.pi))*rt_s*Gb
	return T4s

##### Dielectric losses

def Cable_C(rp,Di, Ds):
	### cp : relative permittivity of the insulation
	### Di : external diameter of insulation (exclude screen)
	### Ds : OD of conductor screen	
	B=18*ln(Ds/Di)
	C=abs((rp/B)*(10**-9))
	return C

def Wd(Cable_C, f, U0, tan_0):
	### Cable_C: capacitance for circular conductor
	### U0 : voltage to earth
	### tan_0 : see IEC table 3
	w=2*math.pi*f
	Wd=w*Cable_C*(U0**2)*tan_0
	return Wd

##### Loss factor for sheath and screen
def lm1_1_flat_trans(Rs, R_AC, f, sheath_ID, sheath_OD, s1):
	w=2*math.pi*f
	a=2*(2**(1/3))
	d=(sheath_ID+sheath_OD)/2
	### per unit length sheath reactance
	X1=(2*w*ln(a*(s1/d)))*(10**-7)
	A=Rs/R_AC
	lm1_1=A*(1/(1+(Rs/X1)**2))
	return lm1_1

def lm1_1_trefoil(Rs, R_AC, f, sheath_ID, sheath_OD, s1):
	w=2*math.pi*f
	d=(sheath_ID+sheath_OD)/2
	### per unit length sheath reactance
	X=(2*w*ln(2*(s1/d)))*(10**-7)
	A=Rs/R_AC
	lm1_1=A*(1/(1+(Rs/X)**2))
	return lm1_1


### Eddy current losses factor

def B1(f, Rss_20, t_s, F_a,):
	w=2*math.pi*f
	Rss_t=Rss_20*(1+F_a*(t_s-20))
	B1= ((4*math.pi*w)/(Rss_t*10**7))**0.5
	return B1

def gs(B1,sheath_ID, sheath_OD):
	thick=(sheath_OD-sheath_ID)/2
	A=(thick/sheath_OD)**1.74
	B=B1*sheath_OD/1000-1.6
	gs=1+A*B
	return gs

def m(f,Rs):
	w=2*math.pi*f
	m=(w/Rs)/10**7
	return m

def lm0_trefoil(m, sheath_ID, sheath_OD, overall_OD):
	A=(m**2)/(1+m**2)
	d=(sheath_ID+sheath_OD)/2
	s=overall_OD
	B=(d/(2*s))**2
	lm_0=3*A*B
	return lm_0

def lm0_flat(m, sheath_ID, sheath_OD, s1):
	A=(m**2)/(1+m**2)
	d=(sheath_ID+sheath_OD)/2
	s=s1
	B=(d/(2*s))**2
	lm_0=6*A*B
	return lm_0

def delta_1_trefoil(m, sheath_ID, sheath_OD, overall_OD):
	if m<=0.1:
		delta_1=0
	else:
		d=(sheath_ID+sheath_OD)/2
		s=overall_OD
		A=1.14*(m**2.45)+0.33
		B=(d/(2*s))**(0.92*m+1.66)
		delta_1=A*B
	return delta_1

def delta_1_flat(m, sheath_ID, sheath_OD, s1):
	if m<=0.1:
		delta_1=0
	else:
		d=(sheath_ID+sheath_OD)/2
		s=s1
		A=0.86*(m**3.08)
		B=(d/(2*s))**(1.4*m+0.7)
		delta_1=A*B
	return delta_1

def lm2 (R_AC,Rs,gs,lm0,delta_1,delta_2,B1,t_s):
	A=Rs/R_AC
	B=gs*lm0*(1+delta_1+delta_2)
	C=((B1*t_s)**4)
	D=12*1e12
	lm2=A*(B+C/D)
	return lm2/10000


